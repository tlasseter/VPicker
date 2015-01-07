
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

public class StsDefineGrid extends StsWizardStep
{
    StsModelGridPanel panel;
    StsHeaderPanel header;
    StsSurface[] surfaces;
    StsGridDefinition defaultGridDef = null;
    byte boundingBoxType = UNION_BOX;

    static final byte UNION_BOX = 1;  // a bounding box around all the grids (union)
    static final byte INTERSECTION_BOX = 2; // bounding box shared by all grids (intersection)
    static final byte SPECIFIED_BOX = 3; // bounding box defined by one of the surfaces

    private StsMatrix4f xyMatrix = null; // transform world coordinates XY to local rotated coordinates
    private StsMatrix4f inverseXYMatrix = null; // inverse transform local rotated XY coordinates to world coordiantes

    public StsDefineGrid(StsWizard wizard)
    {
    	super(wizard, new StsModelGridPanel(), null, new StsHeaderPanel());
        panel = (StsModelGridPanel) getContainer();
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Horizon Construction");
        header.setSubtitle("Define Grid");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");                                        
//        panel.setPreferredSize(new Dimension(400,300));
        header.setInfoText(wizardDialog,"(1) Define a common grid to convert all selected surfaces to.\n" +
                                  "   ***** New grid must be subset of selected horizons. *****\n " +
                                  "   ***** Accepting defaults defines a grid that encompasses all surfaces\n\n" +
                                  "   ***** Minimum X and Y of new grid  ***** \n" +
                                  "   ***** Maximum X and Y of new grid  ***** \n" +
                                  "   ***** Increment in X and Y of new grid ***** \n" +
                                  "   ***** Number of grid cells in X and Y ***** \n" +
                                  "(2) Press the Next>> Button.");
    }

	public void setSurfaces(StsSurface[] surfaces)
    {
    	this.surfaces = surfaces;
        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getStsClass(StsSurface.class);
        defaultGridDef = surfaceClass.checkComputeUnionGrid(surfaces);
        panel.setGridDefinition(defaultGridDef);
    }

    public void setModel(StsActionManager actionManager) { this.model = model; }

    public boolean start()
    {
		surfaces = ((StsHorizonsWizard)wizard).getSurfaces();
		StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getStsClass(StsSurface.class);
        defaultGridDef = surfaceClass.checkComputeUnionGrid(surfaces);
		panel.setGridDefinition(defaultGridDef);
        return true;
    }

    public boolean end()
    {
        return true;
    }

    private void createTransformMatrices(float xOrigin, float yOrigin, float angle)
    {
        // make matrix to transform from world to index cube coordinates

        xyMatrix = new StsMatrix4f();
        inverseXYMatrix.setIdentity();
        xyMatrix.translate(-xOrigin, -yOrigin, 0.0f);
        if(angle != 0.0) xyMatrix.rotZ(-angle);

        // make inverse-matrix to transform from index cube to world coordinates

        inverseXYMatrix = new StsMatrix4f();
        inverseXYMatrix.setIdentity();
        if(angle != 0.0) inverseXYMatrix.rotZ(angle);
        inverseXYMatrix.translate(xOrigin, yOrigin, 0.0f);
    }

    public StsGridDefinition getGridDefinition() { return panel.getGridDefinition(); }
}
