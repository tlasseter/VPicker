
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Export;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.RFUI.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsExportCategoricalFacies extends StsExportZoneLogProperties
{
    static private final String HELP_TARGET = "model.exportFacies";
    static private final String HELP_FILE = "resframe.html";

    private StsCategoricalFacies categoricalFacies = null;
    private int nCategories;
    private boolean okay = true;

 	public StsExportCategoricalFacies(StsActionManager actionManager)
    {
        super(actionManager);
        statusArea.setTitle("Export Categorical Properties:");
        propertyType = CATEGORICAL_TYPE;
    }

    /** override base class query */
    protected String userQuery()
    {
        String logName = askForLogName("Select a categorical facies log curve:",
                "No categorical facies log curve selected.");
        if( logName == null ) return null;

        // if (!dummyTest()) return false;
        if (!buildCategoricalFacies(logName))
        {
            logMessage("Retrieval of log lithology values failed.");
            return null;
        }
        StsCategoricalFaciesDialog d = new StsCategoricalFaciesDialog(null, true,
                categoricalFacies, this);
        d.setLocationRelativeTo(model.win3d);
        d.setVisible(true);
        if (!okay)
        {
            logMessage("Unable to proceed with categorical facies export.");
            return null;
        }
        return logName;
    }


    /* build the categorical facies */
    private boolean buildCategoricalFacies(String logName)
    {
        nCategories = this.getNumberOfCategories();
        if (nCategories == 0) return false;
        try { categoricalFacies = new StsCategoricalFacies(nCategories); }
        catch (StsException e) { return false; }

        // go thru logs and set lithologies
        try
        {
            Iterator wells = model.getObjectIterator(StsWell.class);
            while(wells.hasNext())
            {
                StsWell well = (StsWell)wells.next();
                StsLogCurve logCurve = well.getLastLogCurveOfType(logName);
                if (logCurve != null)
                {
                    float[] values = logCurve.getValuesFloatVector().getValues();
                    categoricalFacies.addLithologies(values);
                }
            }
        }
        catch (NullPointerException e) { return false; }

        return true;
    }

    /* get number of categories */
    private int getNumberOfCategories()
    {
        int nCategories = 4;
        StsSetTextDialog dialog = new StsSetTextDialog(null,
                "Set number of categories (2-9)", (new Integer(nCategories)).toString());
        dialog.setHelpTarget(HELP_TARGET);
        dialog.setHelpFile(HELP_FILE);
        dialog.setLocationRelativeTo(model.win3d);
        dialog.setVisible(true);
        String string = dialog.getText();
        if (string == null)  // cancel
        {
            statusArea.textOnly();
            logMessage("Number of categories unspecified.");
            return 0;
        }
        nCategories = Integer.parseInt(string);
        if (nCategories < 2)
        {
            logMessage("Too few categories specified, setting to 2.");
            nCategories = 2;
        }
        else if (nCategories > 9)
        {
            logMessage("Too many categories specified, setting to 9.");
            nCategories = 9;
        }
        return nCategories;
    }

    protected void writeHeader(String logName) throws IOException
    {
        super.writeHeader(logName);
        out.write("# Number of facies per zone" + NL);
        out.write("GLOBAL" + NL);
        out.write(nCategories + NL);
    }

    /** override base class calculation of layer values */
    protected float[] getSubZoneValues(StsLogCurve logCurve, StsWellZone wellZone)
    {
        float[] values = wellZone.getSubZoneCategoricalAverages(logCurve,
                categoricalFacies);
        if (values == null) return null;

        // fix null values
/*
        for (int i=0; i<values.length; i++)
        {
            if ((int)values[i] == StsCategoricalFacies.INVALID_FACIES ||
                values[i] == StsParameters.nullValue)
            {
                values[i] = StsParameters.HDFnullValue;
            }
        }
*/
        return values;
    }

    public void setOkay(boolean okay) { this.okay = okay; }

    // dummy test
    private boolean dummyTest()
    {
        try { categoricalFacies = new StsCategoricalFacies(5); }
        catch (StsException e)
        {
            logMessage("Couldn't create categorical facies.");
            return false;
        };
        final float[] lithologies = { 1, 2, 3, 4 };
        categoricalFacies.addLithologies(lithologies);

        categoricalFacies.set(1, 1);
        //final int[] facies1 = { 4, 1 };
        //categoricalFacies.addReplace(1, facies1);
        categoricalFacies.set(2, 2);
        categoricalFacies.set(3, 3);
        categoricalFacies.set(4, 4);

        categoricalFacies.print();

        return true;
    }


}
