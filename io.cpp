#include <iostream>
#include <fstream>
#include <vector>

struct Data {
    int id;
    float value;
};

int main() {
    std::ofstream outfile("datafile.bin", std::ios::binary);

    if (!outfile.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return 1;
    }

    // ヘッダの書き込み
    int numElements = 3;  // 例として
    float version = 1.0f;
    outfile.write(reinterpret_cast<char*>(&numElements), sizeof(numElements));
    outfile.write(reinterpret_cast<char*>(&version), sizeof(version));

    // データの書き込み
    std::vector<Data> dataList = {{1, 10.5f}, {2, 20.6f}, {3, 30.7f}};
    for (const Data& data : dataList) {
        outfile.write(reinterpret_cast<char*>(&data), sizeof(data));
    }

    outfile.close();
    return 0;
}

// **********************************************************
#include <iostream>
#include <fstream>
#include <vector>

struct Data {
    int id;
    float value;
};

int main() {
    std::ifstream infile("datafile.bin", std::ios::binary);

    if (!infile.is_open()) {
        std::cerr << "Failed to open file for reading." << std::endl;
        return 1;
    }

    // ヘッダの読み取り
    int numElements;
    float version;
    infile.read(reinterpret_cast<char*>(&numElements), sizeof(numElements));
    infile.read(reinterpret_cast<char*>(&version), sizeof(version));

    std::cout << "Number of elements: " << numElements << std::endl;
    std::cout << "Version: " << version << std::endl;

    // データの読み取り
    std::vector<Data> dataList(numElements);
    for (int i = 0; i < numElements; ++i) {
        infile.read(reinterpret_cast<char*>(&dataList[i]), sizeof(Data));
        std::cout << "Data " << i << ": id=" << dataList[i].id << ", value=" << dataList[i].value << std::endl;
    }

    infile.close();
    return 0;
}


// **********************************************************

#include <iostream>
#include <vector>
#include <tiffio.h>

struct Data {
    int x;
    int y;
    uint8_t value; // 0-255 の値を想定
};

int main() {
    std::vector<Data> dataList = {
        {0, 0, 255},
        {1, 0, 128},
        {0, 1, 64},
        // ... (他のデータ)
    };

    // 画像サイズを定義 (この例では 2x2 としていますが、データに応じて変更してください)
    int width = 2;
    int height = 2;

    uint8_t* image = new uint8_t[width * height];
    for (const Data& data : dataList) {
        image[data.y * width + data.x] = data.value;
    }

    TIFF* out = TIFFOpen("output.tiff", "w");
    if (!out) {
        std::cerr << "Cannot open output.tiff for writing." << std::endl;
        delete[] image;
        return 1;
    }

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

    TIFFWriteEncodedStrip(out, 0, image, width * height);
    TIFFClose(out);
    delete[] image;

    return 0;
}


// **********************************************************

#include <iostream>
#include <vector>
#include <fstream>

struct Data {
    int x;
    int y;
    uint8_t value;
};

int main() {
    std::vector<Data> dataList = {
        {0, 0, 255},
        {1, 0, 128},
        {0, 1, 64},
        {1, 1, 32}
    };

    int width = 2;
    int height = 2;

    // 8ビットグレースケールのバッファを用意
    uint8_t* image = new uint8_t[width * height];
    for (const Data& data : dataList) {
        image[data.y * width + data.x] = data.value;
    }

    std::ofstream tiffFile("output.tif", std::ios::binary);

    // TIFFヘッダの書き込み
    uint16_t endian = 0x4949;  // II (リトルエンディアン)
    tiffFile.write(reinterpret_cast<char*>(&endian), 2);
    uint16_t magic = 42;
    tiffFile.write(reinterpret_cast<char*>(&magic), 2);
    uint32_t ifdOffset = 8;
    tiffFile.write(reinterpret_cast<char*>(&ifdOffset), 4);

    // IFD (Image File Directory) の書き込み
    uint16_t numEntries = 4;
    tiffFile.write(reinterpret_cast<char*>(&numEntries), 2);
    // さまざまなTIFFタグを書き込む...

    // 画像データの書き込み
    tiffFile.write(reinterpret_cast<char*>(image), width * height);

    tiffFile.close();
    delete[] image;

    return 0;
}

