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
    let mut tiff = RasterFile::read("path/to/your/file.tif").unwrap();
    println!("{:?}", tiff.configs);
}
```



