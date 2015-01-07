
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2006
//Author:       TJLasseter
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfacesFromMarkers;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

public class StsBuildSurfaces extends StsWizardStep implements Runnable
{
    public StsProgressPanel panel;
    private StsHeaderPanel header;
    StsMarker[] markers;
    StsSurfacesFromMarkersWizard wizard;

    private boolean isDone = false;
	private boolean canceled = false;

	static private final boolean debugInterpolate = false;

    public StsBuildSurfaces(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsSurfacesFromMarkersWizard)wizard;
        panel = StsProgressPanel.constructorWithCancelButton();
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Create Surfaces from Markers");
        header.setSubtitle("Build the Surfaces");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SurfacesFromMarkers");
        header.setInfoText(wizardDialog,"(1) Once complete, press the Finish Button to dismiss the screen");
    }

    public boolean start()
    {
        run();
        disableNext();
        return true;
    }

    public void run()

    {
        if (canceled)
       {
           success = false;
           return;
       }
        StsModel model = wizard.getModel();
        StsProject project = model.getProject();
//        int nExistingVolumes = model.getObjectList(StsVsp.class).length;
        int n = -1;
        String name = "null";

        try
        {
            StsMarker[] selectedMarkers = wizard.getSelectedMarkers();
			float xInc = wizard.getXInc();
			float yInc = wizard.getYInc();
			StsRotatedGridBoundingBox boundingBox = wizard.boundingBox;
            model.disableDisplay();
            panel.initialize(selectedMarkers.length);
            for(n = 0; n < selectedMarkers.length; n++)
            {
                panel.appendLine("Creating surface from marker: " + selectedMarkers[n].getName());

                // Does it already exist? If so, delete it.
                StsSurface surface = (StsSurface)model.getObjectWithName(StsSurface.class, selectedMarkers[n].getName());
                if(surface != null)
                {
                    panel.appendLine("Deleting existing marker surface: " + selectedMarkers[n].getName());
                    surface.delete();
                }
                // Create surface here..............................
//                surface = StsSurface.constructSurface(selectedMarkers[n].getName(), selectedMarkers[n].getStsColor(), project.getRotatedBoundingBox());
//				surface.setZDomainOriginal(StsProject.TD_DEPTH);
//				surface.initializeSurface();
				StsWellMarker[] wellMarkers = (StsWellMarker[])selectedMarkers[n].getWellMarkers().getCastList();
				surface = leastSqsInterpolateSurface(boundingBox, wellMarkers, selectedMarkers[n].getName(), selectedMarkers[n].getStsColor());
//				surface = modifiedRadiusInterpolateSurface(boundingBox, wellMarkers, selectedMarkers[n].getName(), selectedMarkers[n].getStsColor());
                statusArea.setText("Surface " + surface.getName() + " created successfully.");

                if(surface != null)
                    selectedMarkers[n].setSurface(surface);
                
                panel.setDescription("Surface " + n + " of " + selectedMarkers.length + " created.");
                panel.setValue(n+1);
            }
            isDone = true;
            panel.appendLine("Surface created successfully.");
            panel.setDescription("Complete");
            panel.finished();

        }
        catch(Exception e)
        {
            panel.appendLine("Failed to create surfaces. Error: " + e.getMessage());
            new StsMessage(wizard.frame, StsMessage.WARNING, "Failed to create surfaces");
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            success = false;
            return;
        }
        try
        {
            disableCancel();
            wizard.enableFinish();
            model.enableDisplay();
            model.setActionStatus(StsSurfacesFromMarkersWizard.class.getName(), StsModel.STARTED);
        }
        catch(Exception e)
        {
            panel.appendLine("Failed to create surfaces. Error: " + e.getMessage());
            panel.setDescriptionAndLevel("Exception thrown.", StsProgressBar.ERROR);
            StsException.outputException("StsSurfacesFromMarkers.StsBuildSurfaces.run() failed.", e, StsException.WARNING);
            panel.setDescriptionAndLevel("StsSurfacesFromMarkers.StsBuildSurfaces.run() failed.", StsProgressBar.WARNING);
            success = false;
            return;
        }
    }

	private StsSurface modifiedRadiusInterpolateSurface(StsRotatedGridBoundingBox boundingBox, StsWellMarker[] wellMarkers, String name, StsColor color)
	{
		StsSpiralRadialInterpolation interpolation = new StsSpiralRadialInterpolation(boundingBox, 1, 10);
		interpolation.initialize();
		try
		{
			int nMarkers = wellMarkers.length;
			for(int n = 0; n < nMarkers; n++)
			{
				float[] xyz = wellMarkers[n].getLocation().getXYZorT();
				int nearestRow = boundingBox.getNearestBoundedRowCoor(xyz[1]);
				int nearestCol = boundingBox.getNearestBoundedColCoor(xyz[0]);
                Float zOrT = new Float(xyz[2]);
				interpolation.addDataPoint(zOrT, nearestRow, nearestCol);
			}
			interpolation.run();
			int nRows = interpolation.nRows;
			int nCols = interpolation.nCols;
			float[][] surfaceValues = new float[nRows][nCols];
			for(int row = 0; row < nRows; row++)
			{
				for(int col = 0; col < nCols; col++)
				{
					StsSpiralRadialInterpolation.Weights dataWeights = interpolation.getWeights(row, col);
					if(dataWeights == null || dataWeights.nWeights == 0)
						surfaceValues[row][col] = StsSurface.nullValue;
					else
					{
						int nWeights = dataWeights.nWeights;
						double[] weights = dataWeights.weights;
						Object[] dataObjects = dataWeights.dataObjects;
						double value = 0.0;
						double sumWeight = 0.0;
						for(int n = 0; n < nWeights; n++)
						{
//							StsWellMarker wellMarker = (StsWellMarker)dataObjects[n];
//							float z = wellMarker.getLocation().getZ();
                            float z = ((Float)dataObjects[n]).floatValue();
							if(z != StsParameters.largeFloat)
							{
								value += z*weights[n];
								sumWeight += weights[n];
							}
						}
						if(sumWeight != 0.0)
							surfaceValues[row][col] = (float)(value / sumWeight);
						else
							surfaceValues[row][col] = StsSurface.nullValue;
					}
				}
			}
			return StsSurface.constructSurface(name, color, StsSurface.BUILT, boundingBox, StsProject.TD_TIME, surfaceValues);
		}
		catch(Exception e)
		{
			StsException.outputException("StsBuildSurfacesFromMarkers() failed.", e, StsException.WARNING);
			return null;
		}
	}

	double z, wt, dx, dy;
	double w, wx, wy, wz, wxx, wxy, wyy, wxz, wyz;
    /* ================================================================================================*/
    /* Least Squares Gridding.                                  	                                    */
    /* z = a + b*dx + c*dy                                                                             */
    /* Which means use determinants to solve for constant term a                                        */
    /* We need to solve a system of equations defined as follows                                        */
    /*     wz  =  | w    wx    wy |  a                                                                  */
    /*     wxz =  | wx   wxx  wxy |  b                                                                  */
    /*     wyz =  | wy   wxy  wyy |  c                                                                  */
    /* Using Cramer's rule:                                                                             */
    /*                                                                                                  */
    /* determinant = w*(wxx*wyy - wxy*wxy) + wx*(wy*wxy - wx*wyy) + wy*(wx*wxy - wxx*wy)                */
    /* z0 = a = (wz*(wxx*wyy - wxy*wxy) + wxz*(wxy*wy - wx*wyy) + wyz*(wx*wxy - wxx*wy))/determinant    */
    /* dz/dx = b = w*(wxz*wyy - wyz*wxy) + wz*(wxy*wy - wx*wyy) + wy*(wx*wyz - wxz*wy))/determinant     */
    /* dz/dy = c = w*(wxx*wyz - wxz*wxy) + wx*(wxz*wy - wx*wyz) + wz*(wx*wxy - wxx*wy))/determinant     */
    /*                                                                                                  */
    /* =================================================================================================*/

    private StsSurface leastSqsInterpolateSurface(StsRotatedGridBoundingBox boundingBox, StsWellMarker[] wellMarkers, String name, StsColor color)
	{
		StsSpiralRadialInterpolation interpolation = new StsSpiralRadialInterpolation(boundingBox, 1, 10);
//        StsSurfaceInterpolationWeightedPlane interpolation = new StsSurfaceInterpolationWeightedPlane();
        interpolation.initialize();
		try
		{
			int nMarkers = wellMarkers.length;
			for(int n = 0; n < nMarkers; n++)
			{
				float[] xyz = wellMarkers[n].getLocation().getXYZorT();
				int nearestRow = boundingBox.getNearestBoundedRowCoor(xyz[1]);
				int nearestCol = boundingBox.getNearestBoundedColCoor(xyz[0]);
				interpolation.addDataPoint(wellMarkers[n], nearestRow, nearestCol);
			}
			interpolation.run();
			int nRows = interpolation.nRows;
			int nCols = interpolation.nCols;
			float[][] surfaceValues = new float[nRows][nCols];
			float xMin = boundingBox.xMin;
			float xInc = boundingBox.xInc;
			float yMin = boundingBox.yMin;
			float yInc = boundingBox.yInc;
			float y = yMin;
			for(int row = 0; row < nRows; row++, y += yInc)
			{
				float x = xMin;
				for(int col = 0; col < nCols; col++, x += xInc)
				{
					StsSpiralRadialInterpolation.Weights dataWeights = interpolation.getWeights(row, col);
					if(dataWeights == null || dataWeights.nWeights == 0)
						surfaceValues[row][col] = StsSurface.nullValue;
					else
					{
						int nWeights = dataWeights.nWeights;
						Object[] dataObjects = dataWeights.dataObjects;
						if(nWeights == 1)
						{
							StsWellMarker wellMarker = (StsWellMarker)dataObjects[0];
                            z = wellMarker.getLocation().getZorT();
//							float[] xyz = wellMarker.getLocation().getXYZorT();
//							z = xyz[2];
							surfaceValues[row][col] = (float)z;
						}
						else
						{
							w = 0.0;
							wx = 0.0;
							wy = 0.0;
							wz = 0.0;
							wxx = 0.0;
							wxy = 0.0;
							wyy = 0.0;
							wxz = 0.0;
							wyz = 0.0;

							double[] weights = dataWeights.weights;
//							double weightMin = weights[maxNPoints - 1];

//							StsWellMarker wellMarker = (StsWellMarker)dataObjects[maxNPoints - 1];
//							float[] xyz = wellMarker.getLocation().getXYZorT();

//							dx = (xyz[0] - x)/xInc;
//							dy = (xyz[1] - y)/yInc;
//							double maxDistance = Math.sqrt(dx*dx + dy*dy);

//							weightMin = 1.0/(dx*dx + dy*dy);
							for(int n = 0; n < nWeights; n++)
//							for(int n = 0; n < maxNPoints - 1; n++)
							{
								StsWellMarker wellMarker = (StsWellMarker)dataObjects[n];
//                                z = wellMarker.getLocation().getZ();
							    float[] xyz = wellMarker.getLocation().getXYZorT();
							    z = xyz[2];
//								wt = weights[n] - weightMin;
								dx = (xyz[0] - x)/xInc;
								dy = (xyz[1] - y)/yInc;
								double distance = Math.sqrt(dx*dx + dy*dy);
								wt = 1/distance;
							/*
								double dwt = (distance - maxDistance);
								if(Math.abs(dwt) < 1.e-03)
									wt = 1/distance;
								else
									wt = dwt * dwt / (distance * maxDistance);
							*/
//								wt = 1.0/(dx*dx + dy*dy) - weightMin;
								w += wt;
								wz += wt * z;
								wx += wt * dx;
								wy += wt * dy;
								wxx += wt * dx * dx;
								wxy += wt * dx * dy;
								wyy += wt * dy * dy;
								wxz += wt * dx * z;
								wyz += wt * dy * z;
							}
							if(computeDeterminantOK(w))
							{
								z = (wz * (wxx * wyy - wxy * wxy) + wxz * (wxy * wy - wx * wyy) + wyz * (wx * wxy - wxx * wy)) / determinant;
								if(debugInterpolate)StsMessageFiles.infoMessage("      determinant OK; using gradient. z: " + z);
							}
							else
							{
								z = wz / w;
								if(debugInterpolate)StsMessageFiles.infoMessage("      determinant too small, not using gradient. z: " + z);
							}
							surfaceValues[row][col] = (float)z;
							if(debugInterpolate)StsMessageFiles.infoMessage("      determinant: " + determinant);
						}
					}
				}
			}
//            return StsSurface.constructSurface(name, color, StsSurface.BUILT, boundingBox, StsParameters.TD_DEPTH, surfaceValues);
			return StsSurface.constructSurface(name, color, StsSurface.BUILT, boundingBox, model.getProject().getZDomain(), surfaceValues);
		}
		catch(Exception e)
		{
			StsException.outputException("StsBuildSurfacesFromMarkers() failed.", e, StsException.WARNING);
			return null;
		}
	}

	double determinant;
	static final double minDeterminant = 1.0e-10;

	boolean computeDeterminantOK(double w)
	{
		if(w == 0.0)
		{
			determinant =  0.0;
			return false;
		}
		else
		{
			determinant = w*(wxx*wyy - wxy*wxy) + wx*(wy*wxy - wx*wyy) + wy*(wx*wxy - wxx*wy);
			return Math.abs(determinant) >= minDeterminant;
		}
    }

    public boolean end()
    {
        return true;
    }
}
