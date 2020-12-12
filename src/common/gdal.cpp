#include <paradem/gdal.h>
#include <paradem/raster.h>

#include <iostream>
#include <memory>

bool WriteGeoTIFF(const char* path, int height, int width, void* pData, GDALDataType type, double* geoTransformArray6Eles, double* min, double* max, double* mean, double* stdDev, double nodatavalue) {
    GDALDataset* poDataset;
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

    GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    char** papszOptions = NULL;
    poDataset = poDriver->Create(path, width, height, 1, type, papszOptions);

    if (geoTransformArray6Eles != NULL)
        poDataset->SetGeoTransform(geoTransformArray6Eles);

    GDALRasterBand* poBand;
    poBand = poDataset->GetRasterBand(1);

    poBand->SetNoDataValue(nodatavalue);

    if (min != NULL && max != NULL && mean != NULL && stdDev != NULL) {
        poBand->SetStatistics(*min, *max, *mean, *stdDev);
    }
    poBand->RasterIO(GF_Write, 0, 0, width, height, pData, width, height, type, 0, 0);

    GDALClose((GDALDatasetH)poDataset);
    return true;
}

// read GeoTIFF (float)
bool readGeoTIFF(const char* path, GDALDataType type, Raster<float>& dem) {
    GDALDataset* poDataset;
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    poDataset = (GDALDataset*)GDALOpen(path, GA_ReadOnly);
    if (poDataset == NULL) {
        std::cout << "fail to read file" << std::endl;
        return false;
    }

    GDALRasterBand* poBand;
    poBand = poDataset->GetRasterBand(1);
    GDALDataType dataType = poBand->GetRasterDataType();
    if (dataType != type) {
        std::cout << "fail to read file" << std::endl;
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }

    dem.geoTransforms = std::make_shared<std::vector<double>>(std::vector<double>(6));
    poDataset->GetGeoTransform(&dem.geoTransforms->at(0));
    dem.NoDataValue = poBand->GetNoDataValue();

    if (!dem.init(poBand->GetYSize(), poBand->GetXSize())) {
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }

    poBand->RasterIO(GF_Read, 0, 0, dem.getWidth(), dem.getHeight(), (void*)&dem, dem.getWidth(), dem.getHeight(), dataType, 0, 0);

    GDALClose((GDALDatasetH)poDataset);
    return true;
}

// read GeoTIFF (double)
bool readGeoTIFF2(const char* path, GDALDataType type, Raster<double>& dem) {
    GDALDataset* poDataset;
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    poDataset = (GDALDataset*)GDALOpen(path, GA_ReadOnly);
    if (poDataset == NULL) {
        std::cout << "fail to read file" << std::endl;
        return false;
    }

    GDALRasterBand* poBand;
    poBand = poDataset->GetRasterBand(1);
    GDALDataType dataType = poBand->GetRasterDataType();
    if (dataType != type) {
        std::cout << "fail to read file" << std::endl;
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }
    dem.geoTransforms = std::make_shared<std::vector<double>>(std::vector<double>(6));
    dem.NoDataValue = poBand->GetNoDataValue();
    poDataset->GetGeoTransform(&dem.geoTransforms->at(0));
    if (!dem.init(poBand->GetYSize(), poBand->GetXSize())) {
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }

    poBand->RasterIO(GF_Read, 0, 0, dem.getWidth(), dem.getHeight(), (void*)&dem, dem.getWidth(), dem.getHeight(), dataType, 0, 0);
    GDALClose((GDALDatasetH)poDataset);
    return true;
}

bool readGeoTIFF3(const char* path, GDALDataType type, Raster<uint8_t>& dem) {
    GDALDataset* poDataset;
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
    poDataset = (GDALDataset*)GDALOpen(path, GA_ReadOnly);
    if (poDataset == NULL) {
        std::cout << "fail to read file" << std::endl;
        return false;
    }
    GDALRasterBand* poBand;
    poBand = poDataset->GetRasterBand(1);
    GDALDataType dataType = poBand->GetRasterDataType();
    if (dataType != type) {
        std::cout << "fail to read file" << std::endl;
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }
    dem.geoTransforms = std::make_shared<std::vector<double>>(std::vector<double>(6));
    dem.NoDataValue = poBand->GetNoDataValue();
    poDataset->GetGeoTransform(&dem.geoTransforms->at(0));
    if (!dem.init(poBand->GetYSize(), poBand->GetXSize())) {
        GDALClose((GDALDatasetH)poDataset);
        return false;
    }
    poBand->RasterIO(GF_Read, 0, 0, dem.getWidth(), dem.getHeight(), (void*)&dem, dem.getWidth(), dem.getHeight(), dataType, 0, 0);
    GDALClose((GDALDatasetH)poDataset);
    dem.NoDataValue = 0;
    return true;
}
