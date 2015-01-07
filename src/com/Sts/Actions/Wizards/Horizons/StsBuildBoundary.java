
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

import com.Sts.Actions.Boundary.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsBuildBoundary extends StsWizardStep
{
    StsRadioButtonSelectPanel panel;
    StsHeaderPanel header;

    StsModelSurface[] surfaces = null;
    StsGridDefinition gridDef;
    StsList gridPoints = new StsList(4);
    StsModelSurface surface;
    Class manualBoundaryClass;

    // convenience copies of flags
/*
    static final int VERTICAL_ROW_MINUS = StsParameters.VERTICAL_ROW_MINUS;
    static final int VERTICAL_COL_MINUS = StsParameters.VERTICAL_COL_MINUS;
    static final int VERTICAL_ROW_PLUS = StsParameters.VERTICAL_ROW_PLUS;
    static final int VERTICAL_COL_PLUS = StsParameters.VERTICAL_COL_PLUS;
*/
    static final int NONE = StsParameters.NONE;


    public StsBuildBoundary(StsWizard wizard)
    {
        super(wizard,
        	new StsRadioButtonSelectPanel("Select Boundary Construction Option:",
                new String[] { "Auto build boundary around grid.", "Manual build rectangular boundary.",
                "Manual build polygon boundary." } ), null, new StsHeaderPanel());

        panel = (StsRadioButtonSelectPanel)getContainer();
//        panel.setPreferredSize(new Dimension(400,300));

        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Horizon Construction");
        header.setSubtitle("Build Boundary");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");                                
        header.setInfoText(wizardDialog,"(1) Specify the desired horizon name.\n" +
                                  "(2) Select the surface color (may change from object panel after creation).\n" +
                                  "(3) Enter the Pick Preferences Dialog to adjust graphical refresh rate during picking.\n" +
                                  "   ***** The defaults rate is every pick cycle and should not be adjusted unless ***** \n" +
                                  "   ***** graphics refresh appears to be delayed. *****\n" +
                                  "(4) Press the Next>> Button to proceed to interactive picking.\n");
    }

    public void setSurfaces(StsModelSurface[] surfaces) { this.surfaces = surfaces; }


    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
		int buttonIndexSelected = panel.getButtonIndexSelected();

		switch(buttonIndexSelected)
		{
			case 0:
        	    return constructBoundarySections();
			case 1:
				manualBoundaryClass = StsRectangularBoundary.class;
				return true;
			case 2:
				manualBoundaryClass = StsPolygonBoundary.class;
				return true;
			default:
				return false;
		}
    }

	public boolean constructManualBoundary()
	{
		if(manualBoundaryClass != null)
		    actionManager.startAction(manualBoundaryClass);
		return true;
	}

    private boolean constructBoundarySections()
    {

    	if( surfaces == null ) return false;
    	surface = surfaces[0];
        int ix = surface.getNCols() - 1;
        int iy = surface.getNRows() - 1;
		gridPoints.add(surface.getGridPoint(0,0));
		gridPoints.add(surface.getGridPoint(iy,0));
		gridPoints.add(surface.getGridPoint(iy,ix));
		gridPoints.add(surface.getGridPoint(0,ix));

        /** For each pair of connected points, determine if section is row or col aligned. */
        int i;

        int nLines = 4;
        StsLine[] lines = new StsLine[nLines];
        for(i = 0; i < nLines; i++)
        {
            StsGridPoint gridPoint = (StsGridPoint)gridPoints.getElement(i);
            lines[i] = StsLine.buildVertical(gridPoint, StsParameters.BOUNDARY);
        }

        /** Build sections from vertical pseudos constructed above. */
        int nSides = 4;
        StsSection[] sections = new StsSection[nSides];
		boolean constructionOK = true;
        for(i = 0; i < nSides; i++)
        {
			sections[i] = StsSection.constructor(StsParameters.BOUNDARY, lines[i], lines[(i+1)%nLines]);
		    if(sections[i] == null)
			{
				logMessage("Build Boundary section construction failed for side: " + i +
					"See error log for details.");
				constructionOK = false;
			}
        }
        return constructionOK;
    }
}
