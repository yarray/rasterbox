use super::*;
use utils::Endianness;
use byteorder::{LittleEndian, WriteBytesExt};
use std::f64;
use std::fs::File;
use std::io::prelude::*;
use std::io::Error;
use std::io::ErrorKind;
use std::io::{BufReader, BufWriter, Cursor};
use std::mem;
use std::path::Path;

pub fn read_whitebox(
    file_name: &String,
    configs: &mut RasterConfigs,
    data: &mut Vec<f64>,
) -> Result<(), Error> {
    // read the header file
    // let header_file = file_name.replace(".tas", ".dep");
    let header_file = Path::new(&file_name)
        .with_extension("dep")
        .into_os_string()
        .into_string()
        .unwrap();
    let f = File::open(header_file)?;
    let f = BufReader::new(f);

    for line in f.lines() {
        let line_unwrapped = line.unwrap();
        // println!("{}", line_unwrapped);
        let line_split = line_unwrapped.split(":");
        let vec = line_split.collect::<Vec<&str>>();
        if vec[0].to_lowercase().contains("rows") {
            configs.rows = vec[1].trim().parse::<f32>().unwrap() as usize;
        } else if vec[0].to_lowercase().contains("col") {
            configs.columns = vec[1].trim().parse::<f32>().unwrap() as usize;
        } else if vec[0].to_lowercase().contains("stacks") {
            configs.bands = vec[1].trim().to_string().parse::<u8>().unwrap();
        } else if vec[0].to_lowercase().contains("north") {
            configs.north = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("south") {
            configs.south = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("east") {
            configs.east = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("west") {
            configs.west = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("display min") {
            configs.display_min = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("display max") {
            configs.display_max = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("min")
            && !vec[0].to_lowercase().contains("display")
        {
            configs.minimum = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("max")
            && !vec[0].to_lowercase().contains("display")
        {
            configs.maximum = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("data type") {
            if vec[1].trim().to_lowercase().to_string().contains("double") {
                configs.data_type = DataType::F64;
            } else if vec[1].trim().to_lowercase().to_string().contains("float") {
                configs.data_type = DataType::F32;
            } else if vec[1].trim().to_lowercase().to_string().contains("integer") {
                configs.data_type = DataType::I16;
            } else if vec[1].trim().to_lowercase().to_string().contains("byte") {
                configs.data_type = DataType::U8;
            } else if vec[1].trim().to_lowercase().to_string().contains("i32") {
                configs.data_type = DataType::I32;
            }
        } else if vec[0].to_lowercase().contains("data scale") {
            if vec[1]
                .trim()
                .to_lowercase()
                .to_string()
                .contains("continuous")
            {
                configs.photometric_interp = PhotometricInterpretation::Continuous;
            } else if vec[1]
                .trim()
                .to_lowercase()
                .to_string()
                .contains("categorical")
            {
                configs.photometric_interp = PhotometricInterpretation::Categorical;
            } else if vec[1].trim().to_lowercase().to_string().contains("boolean") {
                configs.photometric_interp = PhotometricInterpretation::Boolean;
            } else if vec[1].trim().to_lowercase().to_string().contains("rgb") {
                configs.photometric_interp = PhotometricInterpretation::RGB;
                configs.data_type = DataType::RGBA32;
            }
        } else if vec[0].to_lowercase().contains("z units") {
            configs.z_units = vec[1].trim().to_string();
        } else if vec[0].to_lowercase().contains("xy units") {
            configs.xy_units = vec[1].trim().to_string();
        } else if vec[0].to_lowercase().contains("projection") {
            configs.projection = vec[1].trim().to_string();
        } else if vec[0].to_lowercase().contains("nodata") {
            configs.nodata = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("preferred palette") {
            configs.palette = vec[1].trim().to_string();
        } else if vec[0].to_lowercase().contains("nonlinearity") {
            configs.palette_nonlinearity = vec[1].trim().to_string().parse::<f64>().unwrap();
        } else if vec[0].to_lowercase().contains("byte order") {
            if vec[1].trim().to_lowercase().contains("little")
                || vec[1].trim().to_lowercase().contains("lsb")
            {
                configs.endian = Endianness::LittleEndian;
            } else {
                configs.endian = Endianness::BigEndian;
            }
        } else if vec[0].to_lowercase().contains("metadata") {
            configs.metadata.push(vec[1].trim().to_string());
        }
    }

    configs.resolution_x = (configs.east - configs.west) / configs.columns as f64;
    configs.resolution_y = (configs.north - configs.south) / configs.rows as f64;

    // read the data file
    // let data_file = file_name.replace(".dep", ".tas");
    let data_file = Path::new(&file_name)
        .with_extension("tas")
        .into_os_string()
        .into_string()
        .unwrap();
    let mut f = File::open(data_file.clone())?;
    //let br = BufReader::new(f);
    // let metadata = fs::metadata(data_file.clone())?;
    // let file_size: usize = metadata.len() as usize;
    // let mut buffer = vec![0; file_size];

    let data_size = if configs.data_type == DataType::F64 {
        8
    } else if configs.data_type == DataType::F32
        || configs.data_type == DataType::I32
        || configs.data_type == DataType::RGBA32
    {
        4
    } else if configs.data_type == DataType::I16 {
        2
    } else {
        // DataType::Byte
        1
    };

    data.reserve(configs.rows * configs.columns);

    let num_cells = configs.rows * configs.columns;
    let buf_size = if num_cells > 10_000_000usize {
        10_000_000usize
    } else {
        num_cells
    };
    while data.len() < num_cells {
        // let mut buffer = vec![0u8; buf_size * data_size];
        let mut buffer = vec![];
        buffer.reserve_exact(buf_size * data_size);
        unsafe {
            buffer.set_len(buf_size * data_size);
        }

        f.read(&mut buffer)?;

        let num_values = if data.len() + buf_size <= num_cells {
            buf_size
        } else {
            num_cells - data.len() + 1
        };

        // read the file's bytes into a buffer
        let mut bor = ByteOrderReader::<Cursor<Vec<u8>>>::new(Cursor::new(buffer), configs.endian);

        match configs.data_type {
            DataType::F64 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_f64()? as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            DataType::F32 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_f32()? as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            DataType::I32 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_i32()? as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            DataType::I16 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_i16()? as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            DataType::U8 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_u8()? as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            DataType::RGBA32 => {
                bor.seek(0);
                for _ in 0..num_values {
                    data.push(bor.read_f32()? as i32 as u32 as f64);
                }
                if data.len() == num_cells {
                    break;
                }
            }
            _ => {
                return Err(Error::new(
                    ErrorKind::NotFound,
                    "Raster data type is unknown.",
                ));
            }
        }
    }

    Ok(())
}

pub fn write_whitebox<'a>(r: &'a mut Raster) -> Result<(), Error> {
    // figure out the minimum and maximum values
    for val in &r.data {
        let v = *val;
        if v != r.configs.nodata {
            if v < r.configs.minimum {
                r.configs.minimum = v;
            }
            if v > r.configs.maximum {
                r.configs.maximum = v;
            }
        }
    }

    if r.configs.display_min == f64::INFINITY {
        r.configs.display_min = r.configs.minimum;
    }
    if r.configs.display_max == f64::NEG_INFINITY {
        r.configs.display_max = r.configs.maximum;
    }

    // Delete the wstat file if it exists
    // let wstat_string = r.file_name.replace(".tas", ".wstat").replace(".dep", ".wstat");
    let wstat_string = Path::new(&r.file_name)
        .with_extension("wstat")
        .into_os_string()
        .into_string()
        .unwrap();
    let wstat_path = Path::new(&wstat_string);
    if wstat_path.exists() {
        match std::fs::remove_file(&wstat_path) {
            Ok(_) => {}  // do nothing
            Err(_) => {} // do nothing
        }
    }

    // Save the header file
    // let header_file = r.file_name.replace(".tas", ".dep");
    let header_file = Path::new(&r.file_name)
        .with_extension("dep")
        .into_os_string()
        .into_string()
        .unwrap();
    let f = File::create(header_file)?;
    let mut writer = BufWriter::new(f);

    let s = format!("Min:\t{}\n", r.configs.minimum);
    writer.write_all(s.as_bytes())?; //.expect("Unable to write data)

    let s = format!("Max:\t{}\n", r.configs.maximum);
    writer.write_all(s.as_bytes())?;

    let s = format!("North:\t{}\n", r.configs.north);
    writer.write_all(s.as_bytes())?;

    let s = format!("South:\t{}\n", r.configs.south);
    writer.write_all(s.as_bytes())?;

    let s = format!("East:\t{}\n", r.configs.east);
    writer.write_all(s.as_bytes())?;

    let s = format!("West:\t{}\n", r.configs.west);
    writer.write_all(s.as_bytes())?;

    let s = format!("Cols:\t{}\n", r.configs.columns);
    writer.write_all(s.as_bytes())?;

    let s = format!("Rows:\t{}\n", r.configs.rows);
    writer.write_all(s.as_bytes())?;

    let s = format!("Stacks:\t{}\n", r.configs.bands);
    writer.write_all(s.as_bytes())?;

    // if r.configs.photometric_interp == PhotometricInterpretation::RGB {
    //     r.configs.data_type = DataType::I32;
    // }

    match r.configs.data_type {
        DataType::F64 | DataType::U32 => {
            if r.configs.photometric_interp != PhotometricInterpretation::RGB {
                // Java doesn't have an unsigned 32-bit integer, so Whitebox only has an I32.
                writer.write_all("Data Type:\tDOUBLE\n".as_bytes())?;
            } else {
                writer.write_all("Data Type:\tI32\n".as_bytes())?;
            }
        }
        DataType::F32 | DataType::RGBA32 | DataType::U16 => {
            writer.write_all("Data Type:\tFLOAT\n".as_bytes())?;
        }
        DataType::I32 => {
            // writer.write_all("Data Type:\tI32\n".as_bytes())?;
            writer.write_all("Data Type:\tFLOAT\n".as_bytes())?;
        }
        DataType::I16 => {
            writer.write_all("Data Type:\tINTEGER\n".as_bytes())?;
        }
        DataType::U8 | DataType::I8 => {
            writer.write_all("Data Type:\tBYTE\n".as_bytes())?;
        }
        _ => {
            return Err(Error::new(
                ErrorKind::NotFound,
                format!(
                    "Raster data type {:?} not supported in this format.",
                    r.configs.data_type
                ),
            ));
        }
    }

    let s = format!("Z Units:\t{}\n", r.configs.z_units);
    writer.write_all(s.as_bytes())?;

    let s = format!("XY Units:\t{}\n", r.configs.xy_units);
    writer.write_all(s.as_bytes())?;

    let s = format!("Projection:\t{}\n", r.configs.projection);
    writer.write_all(s.as_bytes())?;

    match r.configs.photometric_interp {
        PhotometricInterpretation::Continuous => {
            writer.write_all("Data Scale:\tcontinuous\n".as_bytes())?;
        }
        PhotometricInterpretation::Categorical | PhotometricInterpretation::Paletted => {
            writer.write_all("Data Scale:\tcategorical\n".as_bytes())?;
        }
        PhotometricInterpretation::Boolean => {
            writer.write_all("Data Scale:\tBoolean\n".as_bytes())?;
        }
        PhotometricInterpretation::RGB => {
            writer.write_all("Data Scale:\trgb\n".as_bytes())?;
        }
        PhotometricInterpretation::Unknown => {
            writer.write_all("Data Scale:\tcontinuous\n".as_bytes())?;
        }
    }

    let s = format!("Display Min:\t{}\n", r.configs.display_min);
    writer.write_all(s.as_bytes())?;

    let s = format!("Display Max:\t{}\n", r.configs.display_max);
    writer.write_all(s.as_bytes())?;

    if r.configs.palette == String::from("not specified") {
        r.configs.palette = "grey.plt".to_string();
    }
    let s = format!("Preferred Palette:\t{}\n", r.configs.palette);
    writer.write_all(s.as_bytes())?;

    let s = format!("NoData:\t{}\n", r.configs.nodata);
    writer.write_all(s.as_bytes())?;

    if r.configs.endian == Endianness::LittleEndian {
        writer.write_all("Byte Order:\tLITTLE_ENDIAN\n".as_bytes())?;
    } else {
        writer.write_all("Byte Order:\tBIG_ENDIAN\n".as_bytes())?;
    }

    if r.configs.palette_nonlinearity < 0.0 {
        r.configs.palette_nonlinearity = 1.0;
    }
    let s = format!(
        "Palette Nonlinearity:\t{}\n",
        r.configs.palette_nonlinearity
    );
    writer.write_all(s.as_bytes())?;

    for md in &r.configs.metadata {
        let s = format!("Metadata Entry:\t{}\n", md.replace(":", ";"));
        writer.write_all(s.as_bytes())?;
    }

    let _ = writer.flush();

    // write the data file
    // let data_file = r.file_name.replace(".dep", ".tas");
    let data_file = Path::new(&r.file_name)
        .with_extension("tas")
        .into_os_string()
        .into_string()
        .unwrap();
    let f = File::create(&data_file)?;
    let mut writer = BufWriter::new(f);

    // let mut u16_bytes: [u8; 2];
    let mut u32_bytes: [u8; 4];
    let mut u64_bytes: [u8; 8];

    let num_cells: usize = r.configs.rows * r.configs.columns;
    match r.configs.data_type {
        DataType::F64 | DataType::U32 => {
            if r.configs.photometric_interp != PhotometricInterpretation::RGB {
                for i in 0..num_cells {
                    u64_bytes = unsafe { mem::transmute(r.data[i]) };
                    writer.write(&u64_bytes)?;
                }
            } else {
                for i in 0..num_cells {
                    u32_bytes = unsafe { mem::transmute(r.data[i] as u32) };
                    writer.write(&u32_bytes)?;
                }
            }
        }
        DataType::F32 | DataType::U16 => {
            for i in 0..num_cells {
                writer.write_f32::<LittleEndian>(r.data[i] as f32)?;
            }
        }
        DataType::I32 => {
            for i in 0..num_cells {
                writer.write_f32::<LittleEndian>(r.data[i] as f32)?;
            }
        }
        DataType::RGBA32 => {
            for i in 0..num_cells {
                u32_bytes = unsafe { mem::transmute(r.data[i] as u32 as i32 as f32) };
                writer.write(&u32_bytes)?;
            }
        }
        DataType::RGB24 => {
            // The Whitebox raster format doesn't really support a 24-bit RGB;
            // instead use a 32-bit RGBa with saturated alpha channel.
            let mut val: u32;
            let alpha_mask = (255 << 24) as u32;
            for i in 0..num_cells {
                val = alpha_mask | (r.data[i] as u32);
                u32_bytes = unsafe { mem::transmute(val) };
                writer.write(&u32_bytes)?;
            }
        }
        DataType::I16 => {
            for i in 0..num_cells {
                // u16_bytes = unsafe { mem::transmute(r.data[i] as u16) };
                // writer.write(&u16_bytes)?;
                writer.write_i16::<LittleEndian>(r.data[i] as i16)?;
            }
        }
        DataType::U8 | DataType::I8 => {
            for i in 0..num_cells {
                writer.write(&[r.data[i] as u8])?;
            }
        }
        _ => {
            return Err(Error::new(
                ErrorKind::NotFound,
                "Raster data type is unknown.",
            ));
        }
    }

    let _ = writer.flush();

    Ok(())
}
