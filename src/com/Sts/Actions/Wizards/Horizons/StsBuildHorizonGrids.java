
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
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsBuildHorizonGrids extends StsWizardStep implements Runnable
{
    StsStatusPanel status;
    StsHeaderPanel header;
    StsSurface[] surfaces = null;
    StsModelSurface[] horizons = null;
    StsGridDefinition gridDef;

    public StsBuildHorizonGrids(StsWizard wizard)
    {
    	super(wizard, new StsStatusPanel(), null, new StsHeaderPanel());
        status = (StsStatusPanel) getContainer();
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Horizon Construction");
        header.setSubtitle("Build Horizon Grids");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Horizons");                                        
        header.setInfoText(wizardDialog,"(1) Press the Next>> Button when construction is complete.");
        status.setTitle("Building the horizon grids:");
    }

	public void setGridDef(StsGridDefinition def) { this.gridDef = def; }
	public StsModelSurface[] getHorizons() { return horizons; }

	public boolean start()
    {
		run();
		return true;
    }

	/** Create new surfaces of the same size, copying points
	   *  from the old surfaces. Turn off visibility of current surfaces.
	   */
    public void run()
    {
        disableNext();
  		// status.setTitle("Building the horizon grids:");
		surfaces = ((StsHorizonsWizard)wizard).getSurfaces();
		if(surfaces == null) return;
		int nSurfaces = surfaces.length;
		if(nSurfaces == 0) return;
		status.setMinimum(0);
		status.setMaximum(nSurfaces);
		status.setProgress(0);
		status.setText("Initializing ...");

        nSurfaces = surfaces.length;
		horizons = new StsModelSurface[nSurfaces];
		StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
		int nHorizons = 0;
        for(int i = 0; i < nSurfaces; i++)
        {
            StsSurface surface = surfaces[i];
            try
            {
				// get the model surface built from this originalSurface if it exists
				StsModelSurface modelSurface = modelSurfaceClass.getModelSurfaceForOriginal(surface);
	            // new surface
                if(modelSurface == null)
		            status.setText("Creating horizon surface from: " + surface.getName());
                else
                    modelSurface.delete();

                modelSurface = StsModelSurface.constructModelSurfaceFromSurface(surface, modelSurfaceClass, gridDef);
				horizons[nHorizons++] = modelSurface;
//                if(nHorizons > 1) modelSurface.truncateBySurfaceAbove(horizons[nHorizons-2]);
            }
            catch(Exception e)
            {
                StsException.outputException("buildHorizonSurfaces failed.",
                        e, StsException.WARNING);
            }
            surface.setIsVisible(false);
            status.setProgress(i+1); // update progress
	        status.sleep(500);
        }

	    if(nHorizons != nSurfaces)
			horizons = (StsModelSurface[])StsMath.trimArray(horizons, nHorizons);
        // unload IMPORTED surfaces: we are finished, if user redisplays, they can
        // be individually reloaded

        if(nHorizons > 1)
        {
            int firstFullGridIndex = 0;
            for(int n = 0; n < nHorizons; n++)
            {
                if(horizons[n].getOriginalSurface().xyGridSameAs(gridDef))
                {
                    firstFullGridIndex = n;
                    break;
                }
            }
            StsSurface[] thicknessGrids = new StsSurface[nHorizons-1];
            byte zDomain = model.getProject().getZDomain();
            for(int n = 1; n < nHorizons; n++)
               thicknessGrids[n-1] = StsSurface.constructThicknessGrid(horizons[n-1], horizons[n], gridDef, zDomain);

            adjustGridsUsingThickness(horizons, firstFullGridIndex, thicknessGrids);

            for(int n = 0; n < nHorizons; n++)
            {
                horizons[n].saveSurface();
            }
        }

        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsSurface.class);
        surfaceClass.unloadSurfaces();
        success = true;
        enableNext();
        status.setText("Horizon surfaces successfully created.");
    }

    private void adjustGridsUsingThickness(StsSurface[] surfaces, int firstFullGridIndex, StsSurface[] thicknessSurfaces)
    {
        int nRows = horizons[0].nRows;
        int nCols = horizons[0].nCols;
        int nSurfaces = surfaces.length;
        boolean[] surfaceNonNull = new boolean[nSurfaces];
        float[] thicknesses = new float[nSurfaces-1];
        for(int row = 0; row < nRows; row++)
        {
            for(int col = 0; col < nCols; col++)
            {
                boolean hasNonNull = false;
                for(int i = 0; i < nSurfaces; i++)
                {
                    byte nullType = surfaces[i].getPointNull(row, col);
                    surfaceNonNull[i] = (nullType != StsParameters.SURF_BOUNDARY && nullType != StsParameters.SURF_GAP);
                    if(surfaceNonNull[i])
                        hasNonNull = true;
                }

//                if(!hasNonNull)
                    surfaceNonNull[firstFullGridIndex] = true;

                for(int i = 0; i < nSurfaces-1; i++)
                   thicknesses[i] =  thicknessSurfaces[i].getPointZ(row, col);

                adjustUsingNonNullSurfaces(surfaces, surfaceNonNull, thicknesses, row, col);
            }
        }
    }

    private void adjustUsingNonNullSurfaces(StsSurface[] surfaces, boolean[] surfaceNotNull, float[] thicknesses, int row, int col)
    {
        int nSurfaces = surfaces.length;
        boolean botNotNull;
        int botNonNullSurfaceIndex;
        int topNonNullSurfaceIndex = -1;
        float topZ = 0.0f, botZ = 0.0f;

        StsSurface botSurface = surfaces[0];
        botNotNull = surfaceNotNull[0];
        if(botNotNull)
        {
            botNonNullSurfaceIndex = 0;
            botZ = botSurface.getPointZ(row, col);
        }
        else
            botNonNullSurfaceIndex = -1;
        for(int n = 1; n < surfaces.length; n++)
        {
            if(botNotNull)
            {
                topNonNullSurfaceIndex = botNonNullSurfaceIndex;
                topZ = botZ;
            }
            botSurface = surfaces[n];
            botNotNull = surfaceNotNull[n];
            if(botNotNull)
            {
                botNonNullSurfaceIndex = n;
                botZ = botSurface.getPointZ(row, col);
            }
            else
                botNonNullSurfaceIndex = -1;

            if(topNonNullSurfaceIndex == -1 && botNonNullSurfaceIndex != -1)
            {
                float sumThickness = 0.0f;
                for(int i = botNonNullSurfaceIndex-1; i > topNonNullSurfaceIndex; i--)
                {
                    sumThickness += thicknesses[i];
                    float z = botZ - sumThickness;
                    surfaces[i].setPointFilled(row, col, z, StsParameters.SURF_BOUNDARY);
                }
                topNonNullSurfaceIndex = botNonNullSurfaceIndex;
            }
            else if(topNonNullSurfaceIndex != -1 && botNonNullSurfaceIndex != -1 && botNonNullSurfaceIndex > topNonNullSurfaceIndex + 1)
            {
                float totalThickness = 0.0f;
                for(int i = topNonNullSurfaceIndex; i < botNonNullSurfaceIndex; i++)
                {
                    thicknesses[i] = thicknesses[i];
                    totalThickness += thicknesses[i];
                }
                float dZ = botZ - topZ;
                float sumThickness = 0.0f;
                for(int i = topNonNullSurfaceIndex + 1; i < botNonNullSurfaceIndex; i++)
                {
                    sumThickness +=  thicknesses[i-1];
                    float f = sumThickness/totalThickness;
                    float z = topZ + f*dZ;
                    surfaces[i].setPointFilled(row, col, z, StsParameters.SURF_BOUNDARY);
                }
                topNonNullSurfaceIndex = botNonNullSurfaceIndex;
            }
            else if(n == nSurfaces-1 && botNonNullSurfaceIndex == -1)
            {
                if(topNonNullSurfaceIndex != -1)
                {
                    float sumThickness = 0.0f;
                    for(int i = topNonNullSurfaceIndex; i < nSurfaces - 1; i++)
                    {
                        sumThickness += thicknesses[i];
                        float z = topZ + sumThickness;
                        surfaces[i + 1].setPointFilled(row, col, z, StsParameters.SURF_BOUNDARY);
                    }
                }

                else
                {
                    float sumThickness = 0.0f;
                    for(int i = 0; i < nSurfaces - 1; i++)
                    {
                        sumThickness += thicknesses[i];
                        float z = topZ + sumThickness;
                        surfaces[i + 1].setPointFilled(row, col, z, StsParameters.SURF_BOUNDARY);
                    }
                }
//                    StsException.systemError(this, "adjustUsingNonNullSurfaces", " row: " + row + " col: " + col +
//                            " for last surface which is null, there is no previous non-null surface from which to apply thicknesses.");
            }
        }
    }

    public boolean end() { return success; }
}
