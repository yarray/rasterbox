use super::byte_order_reader::Endianness;
use byteorder::{BigEndian, LittleEndian, WriteBytesExt};
use std::io::prelude::*;
use std::io::Error;

pub struct ByteOrderWriter<W: Write> {
    is_le: bool,
    writer: W,
    num_bytes_written: usize,
}

impl<W: Write> ByteOrderWriter<W> {
    pub fn new(writer: W, byte_order: Endianness) -> ByteOrderWriter<W> {
        let is_le = byte_order == Endianness::LittleEndian;
        ByteOrderWriter::<W> {
            writer: writer,
            is_le: is_le,
            num_bytes_written: 0,
        }
    }

    pub fn write_bytes(&mut self, bytes: &[u8]) -> Result<(), Error> {
        self.num_bytes_written += bytes.len();
        self.writer.write_all(bytes)
    }

    pub fn write_u16(&mut self, value: u16) -> Result<(), Error> {
        self.num_bytes_written += 2;
        if self.is_le {
            self.writer.write_u16::<LittleEndian>(value)
        } else {
            self.writer.write_u16::<BigEndian>(value)
        }
    }

    pub fn write_u32(&mut self, value: u32) -> Result<(), Error> {
        self.num_bytes_written += 4;
        if self.is_le {
            self.writer.write_u32::<LittleEndian>(value)
        } else {
            self.writer.write_u32::<BigEndian>(value)
        }
    }

    pub fn write_u64(&mut self, value: u64) -> Result<(), Error> {
        self.num_bytes_written += 8;
        if self.is_le {
            self.writer.write_u64::<LittleEndian>(value)
        } else {
            self.writer.write_u64::<BigEndian>(value)
        }
    }

    pub fn write_f64(&mut self, value: f64) -> Result<(), Error> {
        self.num_bytes_written += 8;
        if self.is_le {
            self.writer.write_f64::<LittleEndian>(value)
        } else {
            self.writer.write_f64::<BigEndian>(value)
        }
    }

    /// Returns the number of bytes written
    pub fn len(&mut self) -> usize {
        self.num_bytes_written
        // // self.writer.stream_len().unwrap() as usize
        // self.writer.seek(SeekFrom::End(0)).unwrap() as usize + 1
    }

    pub fn get_inner(&mut self) -> &W {
        &self.writer
    }
}
