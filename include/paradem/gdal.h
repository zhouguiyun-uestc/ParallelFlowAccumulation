#ifndef PARADEM_GDAL_H
#define PARADEM_GDAL_H

#include <paradem/raster.h>
#include <gdal_priv.h>

bool  WriteGeoTIFF(const char* path, int height, int width, void* pData, GDALDataType type, double* geoTransformArray6Eles,
	double* min, double* max, double* mean, double* stdDev, double nodatavalue);

bool readGeoTIFF(const char* path, GDALDataType type, Raster<float>& dem);
bool readGeoTIFF2(const char* path, GDALDataType type, Raster<double>& dem);
bool readGeoTIFF3(const char* path, GDALDataType type, Raster<uint8_t>& dem);


template <class T>
bool readTif(const char* path, Raster<T>& tif) {
	GDALDataset *poDataset;
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

	tif.geoTransforms = std::make_shared<std::vector<double>>(std::vector<double>(6));
	tif.NoDataValue = poBand->GetNoDataValue();
	poDataset->GetGeoTransform(&tif.geoTransforms->at(0));
	if (!tif.init(poBand->GetYSize(), poBand->GetXSize())) {
		GDALClose((GDALDatasetH)poDataset);
		return false;
	}
	poBand->RasterIO(GF_Read, 0, 0, tif.getWidth(), tif.getHeight(), (void*)&tif, tif.getWidth(), tif.getHeight(), dataType, 0, 0);
	GDALClose((GDALDatasetH)poDataset);
	return true;
}
#endif 
