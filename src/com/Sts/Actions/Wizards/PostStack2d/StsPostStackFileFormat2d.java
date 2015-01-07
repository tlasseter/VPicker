package com.Sts.Actions.Wizards.PostStack2d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;

public class StsPostStackFileFormat2d extends StsPostStackFileFormat 
{
    public StsPostStackFileFormat2d(StsSeismicWizard wizard)
	{
		super(wizard);
	}

    protected void setTitlesAndInfo()
    {
        header.setTitle("Poststack 2d SegY Definition");
		header.setSubtitle("Check and Set File Formats");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#PostStack2d");         
		header.setInfoText(wizardDialog,"(1) Navigate to the directory containing the SegY Files using the Dir Button.\n" +
										"    All SegY Files in the selected directory will be placed in the left list.\n" +
										"      **** File names must have a .sgy, .SGY, .segy or .SEGY suffix *****\n" +
										"(2) Select files from the left list and place them in the right using the controls between the lists.\n" +
										"      **** All selected files will be scanned and placed in the table at bottom of screen with stats****\n" +
										"      **** from a random scan for review to ensure files are read correctly and disk space is adequate. ****\n" +
										"      **** If scan fails, adjust parameters on this and next screen before proceeding w/ processing ****\n" +
										"(3) Select a pre-defined Segy Format template if one exists for the selected files\n" +
										"      **** Format templates are saved from previous executions of this wizard.\n" +
										"(4) Set the appropriate sample format (defaults to header value), units, endianness and name.\n" +
										"      **** A single volume name is applied to all selected volumes. ****\n" +
										"      **** All files must be the same format to be processed in one run ****\n" +
										"(5) Set the scan percentage (# of samples of the SegY Files scanned) for verifying correct file reading.\n" +
										"      **** Random scanning of samples is used to determine if file is being read correctly.\n" +
										"(6) Specify whether the selected files are time/depth, the format, header sizes and start Z.\n" +
										"      **** Headers can be overridden if values in header are incorrect, review headers before overriding ****\n" +
										"(7) Once all selected file formats are correctly defined press the Next>> Button\n" +
										"      **** Dont worry if value in table are incorrect, next screen allows trace header mapping.****");
	}
}