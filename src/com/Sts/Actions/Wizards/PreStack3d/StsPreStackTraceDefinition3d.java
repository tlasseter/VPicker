package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.PostStack3d.*;
import com.Sts.Actions.Wizards.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Oct 15, 2007
 * Time: 10:31:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackTraceDefinition3d extends StsPostStackTraceDefinition
{

    public StsPreStackTraceDefinition3d(StsWizard wizard)
    {
        super(wizard);
    }

    protected StsPostStackTraceDefinitionPanel constructPanel()
    {
        return new StsPostStackTraceDefinitionPanel3d(wizard, this);
    }

    public StringBuffer getInfoText()
	{
		StringBuffer sb = new StringBuffer();
		sb.append("      **** Review the trace definition table to ensure trace values are being read correctly. ****\n");
		sb.append("      **** Minimum required fields are the inline and crossline numbers, ensure they are correct. ****\n");
		sb.append("      **** Trace header attributes that will be stored on processing will have grey backgrounds. ****\n");
		sb.append("(1) Set the trace header size if non-standard (240 bytes).\n");
		sb.append("(2) To limit the table to only header attributes that are to be stored press the Selected Attributes box.\n");
		sb.append("(3) Select rows of attributes that require changes and set new format, byte locations and scalar application.\n");
		sb.append("      ***** Scalar is in location 69 and is multiplied if positive, divided if negative per SEGY standard ****\n");
		sb.append("(4) Press Update button to save the edits to the selected trace header attribute definition.\n");
		sb.append("      ***** Any changes will result in a re-scan of the headers to see if attribute is being correctly read. ****\n");
		sb.append("(5) Add additional desired attributes by selecting the row and pressing the Add button.\n");
		sb.append("(6) Remove attributes by changing table to selected attributes, selecting attribute and pressing remove button.\n");
		sb.append("      ***** CHANGES WILL APPLY TO ALL THE FILES IN THE TABLE AT THE BOTTOM ON THE SCREEN. ****\n");
		sb.append("(7) Once the trace header information has been reviewed and values set, press the Next>> Button.\n");
		sb.append("      ***** All values in file table should be correct except coordinates prior to moving to next screen. ****\n");
		return sb;
	}

	public String getTitle()
	{
		return "Pre-stack Trace Header Details";
	}
}
