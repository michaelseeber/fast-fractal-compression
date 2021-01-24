# Fractal Mapping



This library deals with input and output styles of images that are of PNM format with <b>mod 8 height and mod 8 width</b>. Code can be modified to handle images of other dimensions (by filling in non mod 8 dimensions sides with tailing whitespace columns/rows). The PNM format used consists of a text file with .pnm extension where the first line contains: `P2 image_width image_height MaxGrey` followed by `image_width * image_height` lines with single integer (0 <= integer <= MaxGrey) representing pixel value. Pixel values in file are recorded to the equivalent way of reading image pixels from left to right and then down one row (like reading a book).


Sample pnm images can be found in `sample_images` directory.

To convert jpg file to pnm format run the following:

```
convert lll.jpeg RESULT.pgm
g++ -std=c++11 converter.cc -o convert
./convert RESULT.pgm sample.pnm
```
converter script is located under `tools`.




<h5>
Run regularFractal:
</h5>
Compress (creates `map_output_extension` file that contains mapping):


```
cd regularFractal
g++ -std=c++11 -o compress compress.cc ../blockImage.cc ../compareImages.cc
./compress <input_file.pnm> <map_output_extension>
```

Decompress (creates `image_output_name_extension.pnm`):

```
cd regularFractal
g++ -std=c++11 -o decompress decompress.cc ../blockImage.cc ../compareImages.cc
./decompress <map_output_extension> <starting_image_pnm> <image_output_name_extension>
```

<h5>
Run regularFractalWithRotation:
</h5>

Compress (creates `map_output_extension` file that contains mapping):

```
cd regularFractalWithRotation
g++ -std=c++11 -o compress compress.cc ../blockImage.cc ../compareImages.cc
./compress <input_file.pnm> <map_output_extension>
```

Decompress (creates `image_output_name_extension.pnm`):

```
cd regularFractalWithRotation
g++ -std=c++11 -o decompress decompress.cc ../blockImage.cc ../compareImages.cc
./decompress <map_output_extension> <starting_image_pnm> <image_output_name_extension>
```

