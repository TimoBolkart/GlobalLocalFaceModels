#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <string>
class FileHandler
{
public:
	//! Get the file extension of a full file name.
	//! e.g. ...\...\Test.txt returns txt
	//! \param sstrFile			full file name
	//! \return file extension
	static std::string getFileExtension(const std::string& sstrFileName);

	//! Checks if a file name exists.
	//! \param sstrFile			full file name
	//! \return true if file exists
	static bool fileExist(const std::string& sstrFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns Test
	//! \param sstrFile			full file name
	//! \return file name if file exists
	static std::string getFileName(const std::string& sstrFullFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns ...\...
	//! \param sstrFile			full file name
	//! \return file path if file exists
	static std::string getFilePath(const std::string& sstrFileName);
};

#endif