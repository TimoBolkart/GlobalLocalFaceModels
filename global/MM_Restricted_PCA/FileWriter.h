#ifndef FILEWRITER_H
#define FILEWRITER_H

#include "DataContainer.h"

#include <string>
#include <vector>

class FileWriter
{
public:
	//! Save geometry file. Extend this to support other file formats.
	//! \param sstrFileName		full file name of the exported file
	//! \param data				data that should be exported
	//! \return true if successful
	static bool saveFile(const std::string& sstrFileName, const DataContainer& data);

private:

	//! Specific writer for off files.
	//! \param sstrFileName		full file name of the exported file
	//! \param data				data that should be exported
	//! \return true if successful
	static bool writeOff(const std::string& sstrFileName, const DataContainer& data);
};

#endif