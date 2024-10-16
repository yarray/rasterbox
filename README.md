# RasterBox

A simple library for reading and writing a few geographical raster file formats including GeoTiff, ARC ASCII Grid and GRASS ASCII Grid.
The library is extracted and slimmed down from the [whitebox-tools](https://github.com/jblindsay/whitebox-tools) project. Only read/write functionality for the formats mentioned above is kept for the sake of simplicity.

## Usage

Add this to your `Cargo.toml`:

```toml
rasterbox = { git = "https://github.com/yarray/rasterbox" }
```

## Example

```rust
use rasterbox::RasterFile;

fn main() {
    let mut tiff = RasterFile::read_geotiff("path/to/your/file.tif").unwrap();
    println!("{:?}", tiff.configs);
    // manipulate the content
    tiff.write_geotiff("path/to/your/out_file.tif").unwrap();
}
```

## Roadmap

The library is experimental and will be continue if only I found it to have
more potential than other pure-rust libraries such as
[georaster](https://github.com/pka/georaster) or
[geotiff](https://github.com/georust/geotiff).

- [ ] Better documentation
- [ ] Read metadata only
- [ ] Read data with spatial filter
- [ ] Read COG from remote using http range request
- [ ] Handle multiple bands for GeoTiff
- [ ] Expose more functions
- [ ] Add more tests
- [ ] Add more examples