package com.Sts.Utilities;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import java.io.*;

public class StsListFilteredFiles
{
	String inputDir;
	String filter;
	public StsListFilteredFiles(String inputDir, String f, String outputDir)
	{
		this.inputDir = inputDir;
		File dirFile = new File(inputDir);
		File[] files;
		filter = f;

		// It is also possible to filter the list of returned files.
		// This example does not return any files that start with `.'.
		FilenameFilter fileFilter = new FilenameFilter()
		{
			public boolean accept(File dir, String name)
			{
				return name.startsWith(filter);
			}
		};

		boolean success = (new File(outputDir)).mkdir();


		String[] filenames = dirFile.list(fileFilter);
		for (int n = 0; n < filenames.length; n++)
		{
			try
			{
				BufferedReader in = new BufferedReader(new FileReader(inputDir + File.separator + filenames[n]));
				BufferedWriter out = new BufferedWriter(new FileWriter(outputDir + File.separator + filenames[n]));
				out.write("S2SWELLNAME\n");
				int prefixLength = filter.length();
				String name = filenames[n].substring(prefixLength);
				out.write(name + "\n");
				String str;
				while ( (str = in.readLine()) != null)
				{
					out.write(str + "\n");
				}
				in.close();
				out.close();
			}
			catch (IOException e)
			{
			}
		}

	}

	public static void main(String[] args)
	{
		StsListFilteredFiles stsListFilteredFiles = new StsListFilteredFiles("I:/dataS2S/Q", "well-td.", "I:/dataS2S/Q/td");
	}
}
