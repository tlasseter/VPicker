package com.Sts.DBTypes;

import edu.mines.jtk.interp.*;
import edu.mines.jtk.dsp.Sampling;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.Wizards.Velocity.*;
import java.lang.Math;

import com.Sts.Utilities.Interpolation.StsSpiralRadialInterpolation;

import static com.Sts.Actions.Wizards.Velocity.StsEditVelocityPanel.SPLINES;


/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsVelocityGrid extends StsRotatedGridBoundingBox  // , StsDistanceTransformInterpolationFace
{
    StsModelSurface topSurface = null;
    StsModelSurface botSurface = null;
    double timeDatum, depthDatum;
    StsSeismicVolume inputVelocityVolume = null;
	float constantIntervalVelocity;
    transient StsSeismicVelocityModel velocityModel;
    transient double[][] velocities;
//    double[][] errorCorrections = null;
    transient public float velMin = largeFloat, velMax = -largeFloat;
    transient double[][] weights;
    transient StsProgressPanel panel;
//    transient int gridRowMin, gridColMin;

    static final float largeFloat = StsParameters.largeFloat;
    static boolean debugDistTrans = false;
    static int debugDistTransRow, debugDistTransCol;
    static final float nullValue = StsParameters.nullValue;
    static float maxAllowedError = 0.001f;

    public StsVelocityGrid()
    {
    }

    private StsVelocityGrid(StsSeismicVelocityModel velocityModel, StsModelSurface topSurface, StsModelSurface botSurface, StsSeismicVolume inputVelocityVolume, float constantIntervalVelocity, StsProgressPanel panel) throws StsException
    {
        super(false);
        this.velocityModel = velocityModel;
        this.timeDatum = velocityModel.timeDatum;
        this.depthDatum = velocityModel.depthDatum;
        this.topSurface = topSurface;
        this.botSurface = botSurface;
        this.inputVelocityVolume = inputVelocityVolume;
		this.constantIntervalVelocity = constantIntervalVelocity;
        this.panel = panel;
        constructGeometry();
    }

    static public StsVelocityGrid construct(StsSeismicVelocityModel velocityModel, StsModelSurface topSurface, StsModelSurface botSurface, StsSeismicVolume inputVelocityVolume, float constantIntervalVelocity, StsProgressPanel panel)
    {
        if (botSurface == null)
        {
            return null;
        }
        try
        {
            return new StsVelocityGrid(velocityModel, topSurface, botSurface, inputVelocityVolume, constantIntervalVelocity, panel);
        }
        catch (Exception e)
        {
            StsException.outputException("StsVelocityGrid.constructor() failed for botSurface " + botSurface.getName(),
                e, StsException.WARNING);
            return null;
        }
    }

    private StsVelocityGrid(StsSeismicVelocityModel velocityModel, float timeDatum, float depthDatum, StsModelSurface botSurface, StsSeismicVolume inputVelocityVolume, float constantIntervalVelocity, StsProgressPanel panel) throws StsException
    {
        super(false);
        this.velocityModel = velocityModel;
        this.timeDatum = timeDatum;
        this.depthDatum = depthDatum;
        this.botSurface = botSurface;
        this.inputVelocityVolume = inputVelocityVolume;
		this.constantIntervalVelocity = constantIntervalVelocity;
        this.panel = panel;
        constructGeometry();
//		constructIntervalVelocities();
    }

    static public StsVelocityGrid construct(StsSeismicVelocityModel velocityModel, float timeDatum, float depthDatum, StsModelSurface botSurface, StsSeismicVolume inputVelocityVolume, float constantIntervalVelocity, StsProgressPanel panel)
    {
        if (botSurface == null)
        {
            return null;
        }
        try
        {
            return new StsVelocityGrid(velocityModel, timeDatum, depthDatum, botSurface, inputVelocityVolume, constantIntervalVelocity, panel);
        }
        catch (Exception e)
        {
            StsException.outputException(
                "StsVelocityGrid.constructor() failed for bottom time Surface " + botSurface.getName(),
                e, StsException.WARNING);
            return null;
        }
    }

    private void constructGeometry()
    {
        addVolumeBoundingBox(botSurface);
//        gridRowMin = Math.round(velocityVolume.getRowCoor(yMin));
//        gridColMin = Math.round(velocityVolume.getColCoor(xMin));
    }


    public void reconstructIntervalVelocities(StsProgressPanel panel, float markerFactor, int gridType) throws StsException

    {
        this.panel = panel;

        constructIntervalVelocities(markerFactor, gridType);
 //       fieldChanged("velocities", velocities);
    }

    /** If constantIntervalVelocity == 0.0f, then interval velocities must be computed from surfaces and associated markers.
     * @throws StsException if marker for surface is null or the number of wellMarkers for the marker is 0 and
     * the intervalVelocity is 0.0.
     */

    public void constructIntervalVelocities(float markerFactor /* power factor not radius */, int gridType) throws StsException

    {
        StsMarker marker = botSurface.getMarker();
        if (marker == null && constantIntervalVelocity == 0.0f)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, botSurface.getName() + " has no associated well marker and no non-zero interval velocity.");
            constantIntervalVelocity = 1.0f;
        }
        if(constantIntervalVelocity == 0.0f)
        {
            int nWellMarkers = marker.getWellMarkers().getSize();
            if (nWellMarkers == 0)
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, marker.getName() + " for surface " + botSurface.getName() + " has no markers and no non-zero interval velocity.");
                constantIntervalVelocity = 1.0f;
                marker = null;
            }
        }

        weights = new double[nRows][nCols];
        velocities = new double[nRows][nCols];
         /*
                double[] markerErrorCorrections = null;
                if (inputVelocityVolume != null)
                {
                    errorCorrections = new double[nRows][nCols];
                    markerErrorCorrections = new double[nWellMarkers];
                }
         */
        // For each marker, given its xy location, find the corresponding time
        // on the associated timeSurface.  This time corresponds to the depth
        // of the marker.  With the time and depth at the same xy of the surface
        // above, we can compute the interval velocity at that xy.
//        StsGridPoint gridPoint = new StsGridPoint(botSurface);
		if(constantIntervalVelocity == 0.0f)
			computeFromMarkers(marker, markerFactor, gridType);

		else
        {
            if(marker != null)
                adjustMarkersWithConstantIntervalVelocity(marker);
            computeWithConstantIntervalVelocity();
        }
    }


    /**
     * Returns the parameters 'm' and 'b' for an equation y = mx + b, fitted to
     * the data using ordinary least squares regression. The result is returned
     * as a double[], where result[0] --> b, and result[1] --> m.
     *

     *
     * @return The parameters.
     */
    public static double[] getRegression(float[]xvals, float[]yvals)
    {

        int n = xvals.length;
        if (n < 2) {
            throw new IllegalArgumentException("Not enough data.");
        }
        if (yvals.length != n)
        {
            throw new IllegalArgumentException("mismatched data.");
        }

        double sumX = 0;
        double sumY = 0;
        double sumXX = 0;
        double sumXY = 0;
        for (int i = 0; i < n; i++) {
            double x = (double) xvals[i];
            double y = (double) yvals[i];
            sumX += x;
            sumY += y;
            double xx = x * x;
            sumXX += xx;
            double xy = x * y;
            sumXY += xy;
        }
        double sxx = sumXX - (sumX * sumX) / n;
        double sxy = sumXY - (sumX * sumY) / n;
        double xbar = sumX / n;
        double ybar = sumY / n;

        double[] result = new double[2];
        if (sxx == 0)
        {
           result[1]=0;
           result[0] = ybar;
        }
        else
        {
        result[1] = sxy / sxx;
        result[0] = ybar - result[1] * xbar;
        }
        return result;

    }



    // Incomplete implementation of spiral 1/R-sq interpolator
/*
	private StsSurface computeFromMarkersNew(StsMarker marker)
	{
		double botMarkerT, botMarkerZ;
		double topSurfaceT = 0, topSurfaceZ = 0;

		try
		{
            StsWellMarker[] wellMarkers = (StsWellMarker[])marker.getWellMarkers().getCastList();
		    int nWellMarkers = wellMarkers.length;
		    double[] markerIntervalVelocities = new double[nWellMarkers];

            StsSpiralRadialInterpolation interpolation = new StsSpiralRadialInterpolation(this, 1, 10);
		    interpolation.initialize();

			for(int n = 0; n < nWellMarkers; n++)
			{
                StsWellMarker wellMarker = wellMarkers[n];
                if (wellMarker == null) continue;
                StsPoint location = wellMarker.getLocation();
                if(location == null)
                    continue;
                float x = location.getX();
                float y = location.getY();

                float rowF = getRowCoor(y);
                float colF = getColCoor(x);
                if (rowF < 0 || rowF >= nRows)
                {
                    StsMessageFiles.errorMessage("well marker rowF " + rowF + " is out of range of grid.");
                    continue;
                }
                if (colF < 0 || colF >= nCols)
                {
                    StsMessageFiles.errorMessage("well marker colF " + colF + " is out of range of grid.");
                    continue;
                }

                float botSurfaceT = botSurface.interpolateBilinearNoNulls(rowF, colF);
                if(botSurfaceT == botSurface.nullZValue)
                {
                    StsMessageFiles.errorMessage("Failed to find velocity for surface " + botSurface.getName() + " at row " + rowF + " col " + colF);
                    continue;
                }

                // mainDebug check surface rowF and colF
                float surfaceRowF = botSurface.getRowCoor(y);
                float surfaceColF = botSurface.getColCoor(x);
                if (Math.abs(surfaceRowF - rowF) > 0.01f || Math.abs(surfaceColF - colF) > 0.01f)
                {
                    System.err.println("Difference between velocity grid row & col and surface row & col.");
                }

                if(location.hasT())
                {
                    botMarkerT = location.getT();
                    StsMessageFiles.infoMessage("Adjusting well marker " + marker.getName() + " from time " + botMarkerT + " to " + botSurfaceT);
                }
                else
                {
                    botMarkerT = botSurfaceT;
                    location.setT((float)botMarkerT);
                }
                botMarkerZ = location.getZ();

                topSurfaceT = getSurfaceTime(topSurface, rowF, colF);
                topSurfaceZ = getSurfaceDepth(topSurface, rowF, colF);

                double dz = botMarkerZ - topSurfaceZ;
                //TODO need to do a careful units check to determine a minimum acceptable layer thickness
                if(dz < 2.0) continue;
                double dt = botSurfaceT - topSurfaceT;
                if(dt < 1.0)  continue;
                float intervalVelocity = (float)(dz/dt);

                float avgVelocity = (float) ( (botMarkerZ - depthDatum) / (botSurfaceT - timeDatum));
                // don't think it's reasonable to replace interval velocity with average velocity.
                //   Better to bail based on criteria above.  TJL 1/20/10

    //            velocityModel.adjustVelocityRange(intervalVelocity);
                // used for avg velocity volume
                velMin = Math.min(velMin, avgVelocity);
                velMax = Math.max(velMax, avgVelocity);

                markerIntervalVelocities[n] = intervalVelocity;

                System.out.println("Adjusting well marker " + marker.getName() + " well " + wellMarker.well.getName() +
                                   " from time " + botMarkerT + " to " + botSurfaceT);
                System.out.println("    rowF " + rowF + " colF " + colF);

    //            float errorCorrection = nullValue;
                if (inputVelocityVolume != null)
                {
                    float volumeAvgVelocity = (float)velocityModel.scaleMultiplier*inputVelocityVolume.getTrilinearFloatValue(x, y, botSurfaceT);
    //                errorCorrection = avgVelocity / volumeAvgVelocity;
    //                markerErrorCorrections[n] = errorCorrection;
                    System.out.println("    surfaceAvgVelocity " + avgVelocity + " volumeAvgVelocity " +
                                       volumeAvgVelocity); // + " error correction: " + errorCorrection);
                }
                interpolation.addInterpolationPoint(rowF, colF, intervalVelocity); //, errorCorrection);
                //           adjustWell(wellMarker, botSurfaceT, topTimeSurface);
                if (botSurfaceT != botMarkerT)
                {
                    location.setT(botSurfaceT);
                    wellMarker.fieldChanged("location", location);
                }
				int nearestRow = getNearestBoundedRowCoor(y);
				int nearestCol = getNearestBoundedColCoor(x);
                Float zObject = new Float(z);
				interpolation.addDataPoint(zObject, nearestRow, nearestCol);
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
							StsWellMarker wellMarker = (StsWellMarker)dataObjects[n];
							float z = wellMarker.getLocation().getZ();
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
		}
		catch(Exception e)
		{
			StsException.outputException("StsBuildSurfacesFromMarkers() failed.", e, StsException.WARNING);
		}
	}
*/
  private void computeFromMarkers(StsMarker marker, float markerFactor, int gridType)
	{

		double botMarkerT, botMarkerZ;
		double topSurfaceT = 0, topSurfaceZ = 0;

		StsWellMarker[] wellMarkers = (StsWellMarker[])marker.getWellMarkers().getCastList();
		int nWellMarkers = wellMarkers.length;
		double[] markerIntervalVelocities = new double[nWellMarkers];

		int iradius;
        float markerRadius = 3.f;
        float lastIntervalVelocity = 0;

        float[] xvals = new float[nWellMarkers];
        float[] yvals = new float[nWellMarkers];
        float[] zvals = new float[nWellMarkers];
        float[] intvels = new float[nWellMarkers];


        int kvalid=0;
        for (int n = 0; n < nWellMarkers; n++)
            {
			System.out.println("well marker "+n);
            StsWellMarker wellMarker = wellMarkers[n];
            if (wellMarker == null) continue;
            StsPoint location = wellMarker.getLocation();

                if (location == null){
                      StsMessageFiles.errorMessage("wellmarker" + n +" "+wellMarker+" has no location ")  ;
                    continue;
                }

            float x = location.getX();
            float y = location.getY();
            xvals[kvalid] =  x;
            yvals[kvalid] =  y;
            float rowF = getRowCoor(y);
            float colF = getColCoor(x);
			if (rowF < 0 || rowF >= nRows)
			{
				StsMessageFiles.errorMessage("well marker rowF " + rowF + " is out of range of grid.");
				continue;
			}
			if (colF < 0 || colF >= nCols)
			{
				StsMessageFiles.errorMessage("well marker colF " + colF + " is out of range of grid.");
				continue;
			}

            float botSurfaceT = botSurface.interpolateBilinearNoNulls(rowF, colF);
			if(botSurfaceT == botSurface.nullZValue)
			{
				StsMessageFiles.errorMessage("Failed to find velocity for surface " + botSurface.getName() + " at row " + rowF + " col " + colF);
				continue;
			}

            // mainDebug check surface rowF and colF
            float surfaceRowF = botSurface.getRowCoor(y);
            float surfaceColF = botSurface.getColCoor(x);
            if (Math.abs(surfaceRowF - rowF) > 0.01f || Math.abs(surfaceColF - colF) > 0.01f)
            {
                System.err.println("Difference between velocity grid row & col and surface row & col.");
            }

            if(location.hasT())
            {
                botMarkerT = location.getT();
                StsMessageFiles.infoMessage("Adjusting well marker " + marker.getName() + " from time " + botMarkerT + " to " + botSurfaceT);
            }
            else
            {
                botMarkerT = botSurfaceT;
                location.setT((float)botMarkerT);
            }
            botMarkerZ = location.getZ();

            topSurfaceT = getSurfaceTime(topSurface, rowF, colF);
            topSurfaceZ = getSurfaceDepth(topSurface, rowF, colF);

            double dz = botMarkerZ - topSurfaceZ;
            //TODO need to do a careful units check to determine a minimum acceptable layer thickness
           // if(dz < 2.0) continue;
            double dt = botSurfaceT - topSurfaceT;
             //   System.out.println("bot t "+botSurfaceT+" "+topSurfaceT);
            float intervalVelocity = lastIntervalVelocity;
            double mdt = dt < 0 ? -dt : dt;
            if(mdt < 1.0e-6)  {
                //System.out.println("dt =0 "+dt);
                //intvels[n]=0;
                //zvals[n] = (float) botSurfaceT;
                continue;
            };
             //   System.out.println("bot z "+botMarkerZ+" "+topSurfaceZ);
            //if (dz < 0.5f)
            //{
            //       System.out.println("dz < 0.5 "+dz);
            //   continue;
            //}
                if (dz == 0f)
                   continue;
            lastIntervalVelocity = intervalVelocity = (float)(dz/dt);
            intvels[kvalid] = intervalVelocity;
            zvals[kvalid] = (float) botSurfaceT;

            float avgVelocity = (float) ( (botMarkerZ - depthDatum) / (botSurfaceT - timeDatum));
            // don't think it's reasonable to replace interval velocity with average velocity.
            //   Better to bail based on criteria above.  TJL 1/20/10

//            velocityModel.adjustVelocityRange(intervalVelocity);
            // used for avg velocity volume
            velMin = Math.min(velMin, avgVelocity);
            velMax = Math.max(velMax, avgVelocity);

            markerIntervalVelocities[n] = intervalVelocity;


            System.out.println("Adjusting well marker " + marker.getName() + " well " + wellMarker.well.getName() +
                               " from time " + botMarkerT + " to " + botSurfaceT);
            System.out.println("    rowF " + rowF + " colF " + colF);

//            float errorCorrection = nullValue;
            if (inputVelocityVolume != null)
            {
                float volumeAvgVelocity = (float)velocityModel.scaleMultiplier*inputVelocityVolume.getTrilinearFloatValue(x, y, botSurfaceT);
//                errorCorrection = avgVelocity / volumeAvgVelocity;
//                markerErrorCorrections[n] = errorCorrection;
                System.out.println("    surfaceAvgVelocity " + avgVelocity + " volumeAvgVelocity " +
                                   volumeAvgVelocity); // + " error correction: " + errorCorrection);
            }

            if (gridType ==StsEditVelocityPanel.RSQUARE )
             addInterpolationPoint(rowF, colF, intervalVelocity, markerFactor); //, errorCorrection);
            else
                addInterpolationPoint(rowF, colF, intervalVelocity, 2.0f);    // in case other gridders fail
            if (botSurfaceT != botMarkerT)
            {
                location.setT(botSurfaceT);
                wellMarker.fieldChanged("location", location);
            }
            kvalid++;
        }


            float [] kvels = null;
            float [] kxvals = null;
            float [] kyvals = null;
            float [] kzvals = null;

            if (kvalid < 1)
            {
                System.out.println("no valid markers for velocity calc");

            }
            if ((kvalid < 4 && kvalid > 0)  && (gridType != StsEditVelocityPanel.RSQUARE  ))
            {
              normalizeInterpolation();
              System.out.println("degenerate velocity field; < 4 good markers");

            }   else
            {
            kvels=new float[kvalid];
            kxvals=new float[kvalid];
            kyvals=new float[kvalid];
            kzvals=new float[kvalid];

            System.arraycopy(intvels,0,kvels,0,kvalid);
            System.arraycopy(xvals,0,kxvals,0,kvalid);
            System.arraycopy(yvals,0,kyvals,0,kvalid);
            System.arraycopy(zvals,0,kzvals,0,kvalid);

            }

        if (kvalid >= 4 || (gridType == StsEditVelocityPanel.RSQUARE))
        switch (gridType)
        {


            case StsEditVelocityPanel.RSQUARE:
              normalizeInterpolation();
            break;

            case StsEditVelocityPanel.SPLINES:
                double FX = xMin;
                Sampling SX = new Sampling(nCols,xInc,FX);
                double FY = yMin;
                Sampling SY = new Sampling(nRows,yInc,FY);

                SplinesGridder2 spl = new SplinesGridder2(kvels, kxvals, kyvals);
                if (markerFactor > 0.999f) markerFactor=0.999f ;
                spl.setTension((double)markerFactor) ; /* tension [0...1] */
                float [][] vels  = spl.grid(SX, SY);
                if (vels != null)
                for (int i=0; i< nRows; i++)
                  for(int j=0; j< nCols; j++)
                     if (vels[i][j]>=0)
                        velocities[i][j] = vels[i][j];
                     else
                        velocities[i][j]=0;
            break;

            case StsEditVelocityPanel.BLENDED:
                FX = xMin;
                SX = new Sampling(nCols,xInc,FX);
                FY = yMin;
                SY = new Sampling(nRows,yInc,FY);
                BlendedGridder2 bl = new BlendedGridder2(kvels, kxvals, kyvals);
                bl.setSmoothness((double)markerFactor); /* smoothness 0.25 ...n */
                vels =  bl.grid(SX,SY);
                if (vels != null)
                for (int i=0; i< nRows; i++)
                  for(int j=0; j< nCols; j++)
                  if (vels[i][j]>=0)
                     velocities[i][j] = vels[i][j];
                    else
                        velocities[i][j]=0;
                break;

            case StsEditVelocityPanel.v0kz:
                double [] mxb = getRegression(kzvals,kvels);
                double b = mxb[0];
                double m = mxb[1];
                System.out.println("v0 ="+b+" k= "+m);
                for (int i = 0; i < nRows; i++)
                  for (int j=0; j < nCols; j++)
                  {
                     double tv = getSurfaceTime(botSurface, (float) i, (float) j);
                     velocities[i][j]= (float) (b + (m*tv));
                  }
             break;

            case StsEditVelocityPanel.RADIAL:
                            FX = xMin;
                            SX = new Sampling(nCols,xInc,FX);
                            FY = yMin;
                            SY = new Sampling(nRows,yInc,FY);
                            RadialGridder2 bh = new RadialGridder2(new RadialInterpolator2.Biharmonic(), kvels, kxvals, kyvals);
                            //bl.setSmoothness((double)markerFactor); /* smoothness 0.25 ...n */
                            vels =  bh.grid(SX,SY);
                            if (vels != null)
                            for (int i=0; i< nRows; i++)
                              for(int j=0; j< nCols; j++)
                              if (vels[i][j]>=0)
                                 velocities[i][j] = vels[i][j];
                                else
                                    velocities[i][j]=0;
                            break;
            /*
            case StsEditVelocityPanel.RADIAL2:
                                        FX = xMin;
                                        SX = new Sampling(nCols,xInc,FX);
                                        FY = yMin;
                                        SY = new Sampling(nRows,yInc,FY);
                                        RadialGridder2 bh2 = new RadialGridder2(new RadialInterpolator2.WesselBercovici((double)markerFactor),
                                                kvels, kxvals, kyvals);

                                        vels =  bh2.grid(SX,SY);
                                        if (vels != null)
                                        for (int i=0; i< nRows; i++)
                                          for(int j=0; j< nCols; j++)
                                          if (vels[i][j]>=0)
                                             velocities[i][j] = vels[i][j];
                                            else
                                                velocities[i][j]=0;
                                        break;

            */

        }



        weights = null;

        // make minor adjustments to interval velocities so interpolation back
        // to marker location yields correct value of interval velocity at marker
        System.out.println("Error correction checks");
        for (int n = 0; n < wellMarkers.length; n++)
        {
            StsWellMarker wellMarker = wellMarkers[n];
            if (wellMarker == null)
            {
                continue;
            }
            StsPoint location = wellMarker.getLocation();
            float x = location.getX();
            float y = location.getY();
            float rowF = getRowCoor(y);
            float colF = getColCoor(x);
            if(!insideGrid(rowF, colF)) continue;
//            float[] cornerWeights = new float[4];
            int row = StsMath.below(rowF);
			row = (int) rowF ; // jbw
			if(row < 0) row = 0;
            int col = StsMath.below(colF);
			col = (int) colF;
			if(col < 0) col = 0;
            int irow, icol;
			for (irow = row - (int)markerRadius; irow <= row + (int)markerRadius; irow++)
			for (icol = col - (int)markerRadius; icol <= col + (int)markerRadius; icol++)
			{
				double gridIntervalVelocity = StsMath.interpolateGridBilinear(velocities, (float)irow, (float)icol);
				// jbw Tom's method of zero-thickness handling can yield zeros
				if(markerIntervalVelocities[n] <= 0)markerIntervalVelocities[n] = gridIntervalVelocity;
				double adjustment = markerIntervalVelocities[n] - gridIntervalVelocity;
				double dx = (rowF - irow);
				double dy = (colF - icol);
				double dxdy = (dx * dx) + (dy * dy);
				float dxy;
				if (dxdy > 0)
				   dxy = (float)Math.sqrt((double)(dxdy));
			    else
				   dxy = 1;
			    // 1/sqrt(r) influence changed to straight-line interpolator
				// make a circle
				if (dxy > markerRadius) continue;
				// jbw  we should continue until we get a smooth transition
				if (dxy <= 1.f)
					dxy=1.f;
				else
					dxy = 1.f - (float)(dxy)/(markerRadius + 1.f);
				if (dxy > 1.f) dxy=1.f;

				if(!insideGrid(irow, icol)) continue;
				if(StsSeismicVelocityModel.debug)
				   System.out.println("DEBUG: Well marker " + marker.getName() + " well " + wellMarker.well.getName() +
					   "row " + irow + " col " + icol + " adjust by " + adjustment+ "* factor "+dxy);

				velocities[irow][icol] += (adjustment * dxy);

				double adjustedGridVelocity = StsMath.interpolateGridBilinear(velocities, rowF, colF);
				if ((adjustment > .10 * adjustedGridVelocity) || (adjustment < -0.1 * adjustedGridVelocity))
                {
				System.out.println("Large adjustment: "+adjustment+" Well marker " + marker.getName() + " well " + wellMarker.well.getName() +
								   " gridIntervalVelocity " + gridIntervalVelocity +
								   " markerIntervalVelocity " + markerIntervalVelocities[n] +
								   " adjustedGridIntervalVelocity " + adjustedGridVelocity);
				System.out.println("    irow " + irow + " icol " + icol);
                }
			}
        }

        // With these interval velocities, we can now compute the bottom depth surface from
        // the time and depth surfaces above and the bottom time surface.
        float[][] botDepths = botSurface.getAdjPointsZ();
        float[][] botTimes = botSurface.getPointsZ();
        int progressNRows = nRows / 10;

       float minDepth = largeFloat;
       float maxDepth = -largeFloat;
	   int row, col;
       if (topSurface == null)
       {
		   try
		   {
			   for(row = 0; row < nRows; row++)
			   {
				   if(row % progressNRows == 0)
				   {
                       panel.appendLine("Processing row " + row + " of " + nRows);
					   panel.incrementProgress();
				   }
				   for(col = 0; col < nCols; col++)
				   {
					   botDepths[row][col] = (float)(velocities[row][col] *
						   (botTimes[row][col] - timeDatum) + depthDatum);
					   minDepth = Math.min(botDepths[row][col], minDepth);
					   maxDepth = Math.max(botDepths[row][col], maxDepth);
				   }
			   }
		   }
		   catch(Exception e)
		   {
			   StsException.outputException("StsVelocityGrid.computeFromMarkers() failed.", e, StsException.WARNING);
		   }
       }
       else
       {
           float[][] topDepths = topSurface.getAdjPointsZ();
           float[][] topTimes = topSurface.getPointsZ();
		   try
		   {
			   for(row = 0; row < nRows; row++)
			   {
				   if(row % progressNRows == 0)
				   {
                       panel.appendLine("Processing row " + row + " of " + nRows);
					   panel.incrementProgress();
				   }
				   for(col = 0; col < nCols; col++)
				   {
					   botDepths[row][col] = (float)(velocities[row][col] *
						   (botTimes[row][col] - topTimes[row][col]) +
						   topDepths[row][col]);
					   minDepth = Math.min(botDepths[row][col], minDepth);
					   maxDepth = Math.max(botDepths[row][col], maxDepth);
				   }
			   }
		   }
		   catch(Exception e)
		   {
			   StsException.outputWarningException(this, "computeFromMarkers", "Markers on surface: " + marker.getModelSurface().getName(), e);
		   }
	   }
       panel.appendLine("Completed horizon adjustments for markers on surface: " + marker.getModelSurface().getName());
       panel.finished();

       botSurface.setDepthMin(minDepth);
       botSurface.setDepthMax(maxDepth);

       if(StsSeismicVelocityModel.debug) debugCheckMarker(velocityModel, marker);
   }

   /** because we have markers which must be in depth, they will be honored and time will be adjusted */
   private void adjustMarkersWithConstantIntervalVelocity(StsMarker marker)
   {
	   double botMarkerT, botMarkerZ;
	   double topSurfaceT = 0, topSurfaceZ = 0;

	   StsWellMarker[] wellMarkers = (StsWellMarker[])marker.getWellMarkers().getCastList();
	   int nWellMarkers = wellMarkers.length;
	   double[] markerVelocities = new double[nWellMarkers];
	   for (int n = 0; n < nWellMarkers; n++)
	   {
		   StsWellMarker wellMarker = (StsWellMarker) wellMarkers[n];
		   if (wellMarker == null) continue;
		   StsPoint location = wellMarker.getLocation();
		   float x = location.getX();
		   float y = location.getY();
		   float rowF = getRowCoor(y);
		   float colF = getColCoor(x);
		   if (rowF < 0 || rowF >= nRows)
		   {
			   StsMessageFiles.errorMessage("well marker rowF " + rowF + " is out of range of grid.");
			   continue;
		   }
		   if (colF < 0 || colF >= nCols)
		   {
			   StsMessageFiles.errorMessage("well marker colF " + colF + " is out of range of grid.");
			   continue;
		   }

		   float botSurfaceZ = botSurface.interpolateBilinearNoNulls(rowF, colF);
		   if(botSurfaceZ == botSurface.nullZValue)
		   {
			   StsMessageFiles.errorMessage("Failed to find velocity for surface " + botSurface.getName() + " at row " + rowF + " col " + colF);
			   continue;
		   }

		   // mainDebug check surface rowF and colF
		   float surfaceRowF = botSurface.getRowCoor(y);
		   float surfaceColF = botSurface.getColCoor(x);
		   if (Math.abs(surfaceRowF - rowF) > 0.01f || Math.abs(surfaceColF - colF) > 0.01f)
		   {
			   System.err.println("Difference between velocity grid row & col and surface row & col.");
		   }

		   botMarkerT = location.getT();
		   botMarkerZ = location.getZ();
		   topSurfaceT = getSurfaceTime(topSurface, rowF, colF);
		   topSurfaceZ = getSurfaceDepth(topSurface, rowF, colF);
		   double constantVelocityBotMarkerT = (botMarkerZ - topSurfaceZ)/constantIntervalVelocity + topSurfaceT;
           StsMessageFiles.infoMessage("Adjusting well marker " + marker.getName() + " from time " + botMarkerT + " to " + constantVelocityBotMarkerT);
		   System.out.println("Constant interval velocity adjusts well marker " + marker.getName() + " well " + wellMarker.well.getName() +
			   " from time " + botMarkerT + " to " + constantVelocityBotMarkerT);
		   float avgVelocity = (float) ( (botMarkerZ - depthDatum) / (constantVelocityBotMarkerT - timeDatum));
//		   float intervalVelocity = (float) ( (botMarkerZ - topSurfaceZ) / (botSurfaceT - topSurfaceT));

		   // used for avg velocity volume
		   velMin = Math.min(velMin, avgVelocity);
		   velMax = Math.max(velMax, avgVelocity);

		   markerVelocities[n] = constantIntervalVelocity;

//		   System.out.println("Adjusting well marker " + marker.name + " well " + wellMarker.well.getName() +
//							  " from time " + botMarkerT + " to " + botSurfaceT);
		   System.out.println("    rowF " + rowF + " colF " + colF);

//            float errorCorrection = nullValue;
		   if (inputVelocityVolume != null)
		   {
			   float volumeAvgVelocity = inputVelocityVolume.getTrilinearFloatValue(x, y, (float)constantVelocityBotMarkerT);
//                errorCorrection = avgVelocity / volumeAvgVelocity;
//                markerErrorCorrections[n] = errorCorrection;
			   System.out.println("    surfaceAvgVelocity " + avgVelocity + " volumeAvgVelocity " +
								  volumeAvgVelocity); // + " error correction: " + errorCorrection);
		   }
//		   addInterpolationPoint(rowF, colF, constantIntervalVelocity); //, errorCorrection);
		   //           adjustWell(wellMarker, botSurfaceT, topTimeSurface);
		   if (constantVelocityBotMarkerT != botMarkerT)
		   {
               botMarkerT = constantVelocityBotMarkerT;
               location.setT((float)constantVelocityBotMarkerT);
			   wellMarker.fieldChanged("location", location);
		   }
	   }
   }

    /** With this constant interval velocity, we can now compute the bottom depth surface from
	 * the time and depth surfaces above and the bottom time surface.
     */
    private void computeWithConstantIntervalVelocity()
    {
	   float[][] botTimes = botSurface.getPointsZ();
       float[][] botDepths = new float[nRows][nCols];
       botSurface.setAdjPointsZ(botDepths);
       int progressNRows = nRows / 10;

      velMin = constantIntervalVelocity;
      velMax = constantIntervalVelocity;

      float minDepth = largeFloat;
	  float maxDepth = -largeFloat;
	  if (topSurface == null)
	  {
		  for (int row = 0; row < nRows; row++)
		  {
			  if (row % progressNRows == 0)
			  {
				  panel.incrementProgress();
			  }
			  for (int col = 0; col < nCols; col++)
			  {
				  botDepths[row][col] = (float) ((botTimes[row][col] - timeDatum)*constantIntervalVelocity + depthDatum);
                  velocities[row][col] = constantIntervalVelocity;
//				  minTime = Math.min(botTimes[row][col], minTime);
//				  maxTime = Math.max(botTimes[row][col], maxTime);
			  }
		  }
	  }
	  else
	  {
		  float[][] topTimes = topSurface.getPointsZ();
		  float[][] topDepths = topSurface.getAdjPointsZ();
		  for (int row = 0; row < nRows; row++)
		  {
			  if (row % progressNRows == 0)
			  {
				  panel.incrementProgress();
			  }
			  for (int col = 0; col < nCols; col++)
			  {
                  botDepths[row][col] = ((botTimes[row][col] - topTimes[row][col])*constantIntervalVelocity + topDepths[row][col]);
                  velocities[row][col] = constantIntervalVelocity;
			  }
		  }
	  }
  }

    /** Below but not including marker above, adjust the well path to new times
     *  and adjust all markers to the new well path.
     */
    /*
        private void adjustWell(StsWellMarker wellMarker, float newTime,
                                StsModelSurface topTimeSurface)
        {
            StsWellMarker wellMarkerAbove = null;
            StsPoint locationAbove, location;
            StsWell well = wellMarker.getWell();

            if (topTimeSurface != null)
            {
                wellMarkerAbove = well.getMarker(topTimeSurface.getMarkerName());
                if (wellMarkerAbove != null)
                {
                    locationAbove = wellMarkerAbove.getLocation();
                }
                else
                {
                    locationAbove = well.computeGridIntersect(topTimeSurface);
                }
            }
            else
            {
                locationAbove = new StsPoint(5);
                locationAbove.setT(timeMin);
                locationAbove.setZ(depthMin);
                locationAbove.setM(0.0f);
            }
            location = wellMarker.getLocation();
            adjustWellPath(well, newTime, locationAbove, location);
//        adjustWellMarkers(well, locationAbove, location, wellMarkerAbove,
//                          wellMarker);
        }
     */
    /** We have adjusted the time for a marker at the bottom of an interval.
     *  For the well from the top of the interval to the bottom the interval apply
     *  this adjustment using a cubic spline: zero adjustment, zero derivative at the top,
     *  delta-t adjustment, zero derivative at the bottom.
     */
    /*
        private void adjustWellPath(StsWell well, float newTime, StsPoint point0,
                                    StsPoint point1)
        {
            StsObjectRefList lineVertices = well.getLineVertices();
//       StsPoint[] lineVertexPoints = well.getLineVertexPoints();
            float m0 = point0.getM();
            float m1 = point1.getM();
            float t0 = point0.getT();
            float t1 = point1.getT();
            float dt = newTime - t1;
            float dm = m1 - m0;
            float b = 3 * dt / (dm * dm);
            float a = -2 * dt / (dm * dm * dm);
            boolean adjusted = false;

            int nLineVertices = lineVertices.getSize();
            for (int n = 0; n < nLineVertices; n++)
            {
                StsSurfaceVertex lineVertex = (StsSurfaceVertex) lineVertices.getElement(n);
                StsPoint linePoint = lineVertex.getPoint();
                float m = linePoint.getM();
                float t;
                if (m < m0)
                {
                    continue;
                }
                else if (m >= m0 && m <= m1)
                {
                    float dmm = m - m0;
                    float dtt = (b + a * dmm) * (dmm * dmm);
                    if (dtt != 0.0f)
                    {
                        adjusted = true;
                        t = linePoint.getT();
                        linePoint.setT(t + dtt);
                        //                   StsGridSectionPoint surfacePoint = lineVertex.getSurfacePoint();
                        //                   lineVertex.dbFieldChanged("surfacePoint", surfacePoint);
                        //                   surfacePoint.dbFieldChanged("point", linePoint);
                    }
                }
                else if (dt != 0.0f)
                {
                    adjusted = true;
                    t = linePoint.getT();
                    linePoint.setT(t + dt);
//                StsGridSectionPoint surfacePoint = lineVertex.getSurfacePoint();
//                lineVertex.dbFieldChanged("surfacePoint", surfacePoint);
//                   surfacePoint.dbFieldChanged("point", linePoint);
                }
            }
            if (!adjusted)return;
            well.computeXYZPoints();
            well.setInitialized(false); // initialized flag used temporarily until well changes saved to db in StsSeismicVelocityModel
            // on db load, recompute well path because of the above changes
//        currentModel.addMethodCmd(well, "computeXYZPoints", new Object[0], "Velocity edit.");
        }
     */
    /*
        private void adjustWellMarkers(StsWell well, StsPoint locationAbove,
                                       StsPoint markerLocation,
                                       StsWellMarker wellMarkerAbove,
                                       StsWellMarker wellMarker)
        {
            float mdepthAbove = locationAbove.getM();

            StsWellMarker[] wellMarkers = well.getMarkerArray();

            for (int n = 0; n < wellMarkers.length; n++)
            {
                StsPoint location = wellMarkers[n].getLocation();
                float t = location.getT();
                float mdepth = location.getM();
                if (mdepth <= mdepthAbove)
                {
                    continue;
                }
//           float newTime = well.getTimeFromMDepth(mdepth);
                StsPoint newLocation = well.getPointFromMDepth(mdepth);
                float newT = newLocation.getT();
                if (wellMarkers[n] == wellMarkerAbove)
                {
                    System.err.println("Error:  we are adjusting set marker above " +
                                       wellMarkerAbove.getName() + " from " + t + " to " + newT);
                    continue;
                }
                location = newLocation;
            }
        }
     */
    private void addInterpolationPoint(float rowF, float colF, float velocity, float markerFactor) //, float errorCorrection)
    {
        float weight;
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                float distance = (float) Math.sqrt( (row - rowF) * (row - rowF) +
                    (col - colF) * (col - colF));
		// jbw 1/r squared velocity volume interpolation
                if (distance < 1.d)
                {

                    weight = 10.0f;
					//weight = 1.0f;

                }
                else
                {
                    weight = (float) (1.0d / Math.pow( distance,(double) markerFactor));
                    //weight = 1.0f / (float)Math.sqrt((double)distance);
					//weight = 2.0f * distance;
                }
                velocities[row][col] += weight * velocity;
                /*
                    if (inputVelocityVolume != null)
                    {
                        errorCorrections[row][col] += weight * errorCorrection;
                    }
                 */
                weights[row][col] += weight;
				//System.out.println("add vel "+row+" "+col+" "+velocities[row][col]+" "+weight);
            }
        }
    }

    private void normalizeInterpolation()
    {
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                velocities[row][col] /= weights[row][col];
            }
        }
        /*
            if (inputVelocityVolume == null)
            {
                return;
            }

            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    errorCorrections[row][col] /= weights[row][col];
                }
            }
         */
    }

    /*
        public void constructIntervalVelocities() throws StsException
        {
            float botT = 0, botZ = 0, topT = 0, topZ = 0;
            StsMarker marker = botTimeSurface.getMarker();
            if (marker == null)
                throw new StsException(StsException.WARNING, botTimeSurface.getName() + " has no markers.");
            Object[] wellMarkers = marker.getWellMarkers().getArrayList();
            if (wellMarkers == null || wellMarkers.length == 0)
                throw new StsException(StsException.WARNING, botTimeSurface.getName() + " has no markers.");

            distances = new float[nRows][nCols];
            velocities = new float[nRows][nCols];
            for(int row = 0; row < nRows; row++)
                for(int col = 0; col < nCols; col++)
                    distances[row][col] = largeFloat;

            // For each marker, given its xy location, find the corresponding time
            // on the associated timeSurface.  This time corresponds to the depth
            // of the marker.  With the time and depth at the same xy of the surface
            // above, we can compute the interval velocity at that xy.
            StsGridPoint gridPoint = new StsGridPoint(botTimeSurface);
            for (int n = 0; n < wellMarkers.length; n++)
            {
                StsWellMarker wellMarker = (StsWellMarker) wellMarkers[n];
                StsPoint location = wellMarker.getLocation();
                float x = location.getX();
                float y = location.getY();
                float rowF = getRowCoor(y);
                gridPoint.setRowF(rowF);
                float colF = getColCoor(x);
                gridPoint.setColF(colF);
                float botSurfaceT = botTimeSurface.getZ(gridPoint);
                botT = location.getT();
     StsMessageFiles.infoMessage("Adjusting well marker " + marker.name + " from time " + botT + " to " + botSurfaceT);
                botT = botSurfaceT;
                botZ = location.getZ();
                topT = getSurfaceTime(topTimeSurface, gridPoint);
                topZ = getSurfaceDepth(topDepthSurface, gridPoint);
                int row = Math.round(gridPoint.rowF);
                if(row < 0 || row >= nRows)
                {
     System.out.println("well marker row " + row + " is out of range of grid.");
                    continue;
                }
                int col = Math.round(gridPoint.colF);
                if(col < 0 || col >= nCols)
                {
     System.out.println("well marker col " + col + " is out of range of grid.");
                    continue;
                }
                distances[row][col] = 0.0f;
                velocities[row][col] = (botZ - topZ)/(botT - topT);

                location.setT(botT);
            }
            // Now that we have interval velocity points set for each marker (at the nearest row-col),
            // we use a fast 1/r interpolator to fill in an interpolated interval velocity grid.
     StsDistanceTransformInterpolation distanceTransform = new StsDistanceTransformInterpolation(this, false);
            distanceTransform.interpolateDistanceTransform();

            // With these interval velocities, we can now compute the bottom depth surface from
            // the time and depth surfaces above and the bottom time surface.
            float[][] botDepths = botDepthSurface.getPointsZ();
            float[][] botTimes = botTimeSurface.getPointsZ();
            int progressNRows = nRows/10;
            if(topTimeSurface == null)
            {
                for (int row = 0; row < nRows; row++)
                {
                    if(row%progressNRows == 0) panel.incrementProgress();
                    for (int col = 0; col < nCols; col++)
     botDepths[row][col] = velocities[row][col] * (botTimes[row][col] - tMin) + zMin;
                }
            }
            else
            {
                float[][] topDepths = topDepthSurface.getPointsZ();
                float[][] topTimes = topTimeSurface.getPointsZ();
                for(int row = 0; row < nRows; row++)
                {
                    if (row % progressNRows == 0) panel.incrementProgress();
                    for (int col = 0; col < nCols; col++)
                        botDepths[row][col] = velocities[row][col] * (botTimes[row][col] - topTimes[row][col]) +
                            topDepths[row][col];
                }
            }
            debugCheck();
     }
     */
    public void debugCheck(StsSeismicVelocityModel velocityModel)
    {
        StsMarker marker = botSurface.getMarker();
        if(marker != null)
            debugCheckMarker(velocityModel, marker);
        debugCheckSurface(velocityModel);
    }

    private void debugCheckMarker(StsSeismicVelocityModel velocityModel, StsMarker marker)
    {
        double botSurfaceT, botSurfaceZ;
        double botMarkerZ;

        Object[] wellMarkers = marker.getWellMarkers().getArrayList();
        double maxError = 0.0;
        System.out.println("Debug check for interval above surface " + botSurface.getName());
        StsSeismicVolume outputVelocityVolume = velocityModel.velocityVolume;
        float[][] pointsT = botSurface.getPointsZ();
        float[][] pointsZ = botSurface.getAdjPointsZ();
        float velocityScaleMultiplier = (float)velocityModel.scaleMultiplier;
        for (int n = 0; n < wellMarkers.length; n++)
        {
            StsWellMarker wellMarker = (StsWellMarker) wellMarkers[n];
            if (wellMarker == null)
            {
                continue;
            }
            StsPoint location = wellMarker.getLocation();
            float x = location.getX();
            float y = location.getY();
            float rowF = getRowCoor(y);
            float colF = getColCoor(x);
            if(!insideGrid(rowF, colF)) continue;

            //            topSurfaceT = getSurfaceTime(topTimeSurface, gridPoint);
            //            topSurfaceZ = getSurfaceDepth(topDepthSurface, gridPoint);
            botMarkerZ = location.getZ();
            botSurfaceT = botSurface.interpolateBilinearNoNulls(pointsT, rowF, colF);
            botSurfaceZ = botSurface.interpolateBilinearNoNulls(pointsZ, rowF, colF);
            double markerAvgVelocity = (botMarkerZ - depthDatum) / (botSurfaceT - timeDatum);
            double surfaceAvgVelocity = (botSurfaceZ - depthDatum) / (botSurfaceT - timeDatum);
            double error = (surfaceAvgVelocity - markerAvgVelocity) / markerAvgVelocity;
            maxError = Math.max(Math.abs(error), maxError);

            if (error > maxAllowedError || error < -maxAllowedError)
            {
                System.out.println("Avg velocity error at marker " + marker.getName() + " well " + wellMarker.well.getName() +
                                   ". Should be: " +
                                   markerAvgVelocity + " but is " + surfaceAvgVelocity + ". Error: " + error * 100 + "%");
                System.out.println("    rowF " + rowF + " colF " + colF);
            }

            if (outputVelocityVolume != null)
            {
                float outputVolumeAvgVelocity = outputVelocityVolume.getTrilinearFloatValue(x, y, (float) botSurfaceT, false);
//                float inputVolumeAvgVelocity = inputVelocityVolume.getTrilinearFloatValue(x, y, (float) botSurfaceT);
//                 double errorCorrection = surfaceAvgVelocity/inputVolumeAvgVelocity;
//               float correctedVolumeAvgVelocity = (float)(inputVolumeAvgVelocity * errorCorrection);
                double volumeError = (surfaceAvgVelocity - outputVolumeAvgVelocity)/outputVolumeAvgVelocity;
                if (volumeError > maxAllowedError || volumeError < -maxAllowedError)
                {
                    System.out.println("PostStack3d-to-surface avg velocity error at marker " + marker.getName() + " well " +
                                       wellMarker.well.getName() + ". surfaceAvgVelocity " + surfaceAvgVelocity +
                                       " outputVolumeAvgVelocity " + outputVolumeAvgVelocity +
                                       ". Error: " + volumeError * 100 + "%");
                    System.out.println("    rowF " + rowF + " colF " + colF +
                                       " volumeRowF " + outputVelocityVolume.getRowCoor(y) + " volumeColF " + outputVelocityVolume.getColCoor(x));
                }
            }
        }

        System.out.println("Maximum marker interval velocity error: " + 100 * maxError + "%");
        panel.appendLine("Maximum marker interval velocity error: " + 100 * maxError + "%");
    }

    private void debugCheckSurface(StsSeismicVelocityModel velocityModel)
    {
        double topSurfaceT, topSurfaceZ;
        double botSurfaceT, botSurfaceZ;

        double maxError = 0.0;
        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                topSurfaceT = getSurfaceTime(topSurface, row, col);
                topSurfaceZ = getSurfaceDepth(topSurface, row, col);
                botSurfaceT = getSurfaceTime(botSurface, row, col);
                botSurfaceZ = getSurfaceDepth(botSurface, row, col);
                double surfaceIntervalVelocity = (botSurfaceZ - topSurfaceZ) / (botSurfaceT - topSurfaceT);
                double intervalVelocity = velocities[row][col];
                double error = (surfaceIntervalVelocity - intervalVelocity) / intervalVelocity;
                maxError = Math.max(Math.abs(error), maxError);
                if (error > maxAllowedError || error < -maxAllowedError)
                {
                    StsMessageFiles.errorMessage("Interval velocity error at row " + row +
                                                 " col " + col + "." +
                                                 " Interval velocity: " +
                                                 intervalVelocity +
                                                 " surfaces interval velocity " +
                                                 surfaceIntervalVelocity);
                }
            }
        }
        String message = new String("Maximum surface interval velocity error for " + botSurface.getName() + ": " +
                                    100 * maxError + "%");
        StsMessageFiles.infoMessage(message);
        panel.appendLine(message);
    }

    private double getSurfaceTime(StsModelSurface timeSurface, int row, int col)
    {
        if (timeSurface == null)
        {
            return timeDatum;
        }
        return (double) timeSurface.getPointZ(row, col);
    }

    private double getSurfaceDepth(StsModelSurface depthSurface, int row, int col)
    {
        if (depthSurface == null)
        {
            return depthDatum;
        }
        return (double) depthSurface.getAdjPointZ(row, col);
    }

    private double getSurfaceTime(StsModelSurface timeSurface, float rowF, float colF)
    {
        if (timeSurface == null)
        {
            return timeDatum;
        }
        else
        {
            return (double) timeSurface.interpolateTimeNoNulls(rowF, colF);
        }
    }

    private double getSurfaceDepth(StsModelSurface depthSurface, float rowF, float colF)
    {
        if (depthSurface == null)
        {
            return depthDatum;
        }
        else
        {
            return (double) depthSurface.interpolateDepthNoNulls(rowF, colF);
        }
    }

    public int getNRows()
    {
        return nRows;
    }

    public int getNCols()
    {
        return nCols;
    }

    /*
     public float[][] initializeDistances()
     {
      return distances;
     }
     */
    public float[] getDistanceParameters()
    {
        return new float[]
            {
            yInc, xInc};
    }

    public float evaluateDistTransform(int row, int col,
                                       StsDistanceTransformPoint[] points)
    {
        float weight;

        float nearestDistance = largeFloat;
        boolean debug = debugDistanceTransform(row, col);
        float z = 0;
        for (int n = 1; n < 5; n++)
        {

            float distance = points[n].getDistance();
            int nearestRow = points[n].row;
            int nearestCol = points[n].col;
            if (debug)
            {
                System.out.println("    row: " + nearestRow + " col: " +
                                   nearestCol + " distance: " + distance +
                                   " velocity: " +
                                   velocities[nearestRow][nearestCol]);
            }
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                z = (float) velocities[nearestRow][nearestCol];
            }
        }

        float d0 = points[0].distance;

        if (nearestDistance == largeFloat)
        {
            if (debug)
            {
                System.out.println("    no available points to interpolate");
            }
            return d0;
        }

        if (d0 != largeFloat)
        {
            float z0 = (float) velocities[row][col];
            float w0 = 1.0f / d0;
            float w1 = 1.0f / nearestDistance;
            float f = w0 / (w0 + w1);
            if (debug)
            {
                System.out.print("  interpolation between old z: " + z0 +
                                 " and new z: " + z +
                                 " interpolation factor: " + f +
                                 " nearest distance: " + nearestDistance);
            }
            z = f * z0 + (1 - f) * z;
            if (debug)
            {
                System.out.println(" resulting z: " + z);
            }
            if (d0 < nearestDistance)
            {
                nearestDistance = d0;
            }
        }
        else
        {
            if (debug)
            {
                System.out.println("    new z value at this point: " + z +
                                   " nearest distance: " +
                                   nearestDistance);
            }
        }
//        if(wtSum == 0.0f) return largeFloat;

//        float z;
//        if(wtSum == 1.0f)
//            z = zSum;
//        else
//            z = zSum/wtSum;

        velocities[row][col] = z;
        return nearestDistance;
    }

    boolean debugDistanceTransform(int row, int col)
    {
        if (!debugDistTrans)
        {
            return false;
        }
        return row == debugDistTransRow && col == debugDistTransCol;
    }

    public double getVelocity(float rowF, float colF)
    {
        int i = 0, j = 0;
        int n;
        float dx, dy;
        boolean debug = false;
        double z;
        //int nNulls = 0;

        try
        {
            if (velocities == null)
            {
                return nullValue;
            }

            if (!insideGrid(rowF, colF))
            {
                return nullValue;
            }

            i = (int) rowF;
            if (i == nRows - 1)
            {
                i--;
            }
            dy = rowF - i;

            j = (int) colF;
            if (j == nCols - 1)
            {
                j--;
            }
            dx = colF - j;

            if (debug)
            {
                System.out.println("\txInt = " + rowF +
                                   ", dx = " + dx + ", yInt = " + colF +
                                   ", dy = " + dy);
            }

            float weight = 0.0f;
            float zWeighted = 0.0f;
            float w;

            w = (1.0f - dy) * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = (float) velocities[i][j];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i][j] = " + z);
                }
            }

            w = dy * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = (float) velocities[i + 1][j];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i+1][j] = " + z);
                }
            }

            w = (1.0f - dy) * dx;
            if (w > StsParameters.roundOff)
            {
                z = (float) velocities[i][j + 1];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[ijE] = " + z);
                }
            }

            w = dy * dx;
            if (w > StsParameters.roundOff)
            {
                z = (float) velocities[i + 1][j + 1];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i+1][j+1] = " + z);
                }
            }
            return zWeighted / weight;
        }
        catch (Exception e)
        {
            StsException.outputException(
                "StsGrid.interpolateBilnearNulls() failed." +
                " surface: " + getName() + " i: " + i + " j: " + j, e,
                StsException.WARNING);
            return nullValue;
        }
    }

    /** After velocityModel has been computed and on a subsequent DB load,
     *  if interval velocities are being used recompute them from grid time/depth surfaces.
     */

    public boolean initializeIntervalVelocities()
    {
        int row = -1, col = -1;

        try
        {
            float[][] botDepths = botSurface.getAdjPointsZ();
            float[][] botTimes = botSurface.getPointsZ();
            if(velocities == null) velocities = new double[nRows][nCols];

            if (topSurface == null)
            {
                for (row = 0; row < nRows; row++)
                    for (col = 0; col < nCols; col++) {
						velocities[row][col] = (botDepths[row][col] - depthDatum) / (botTimes[row][col] - timeDatum);
						if (velocities[row][col] <= 0)
							System.out.println("neg vel");
					}
            }
            else
            {
                float[][] topDepths = topSurface.getAdjPointsZ();
                float[][] topTimes = topSurface.getPointsZ();
                for (row = 0; row < nRows; row++)
                    for (col = 0; col < nCols; col++) {
                        velocities[row][col] = (botDepths[row][col] - topDepths[row][col]) / (botTimes[row][col] - topTimes[row][col]);
					}

			}

 //           correctIntervalVelocitiesToMarkers();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsVelocityGrid.initializeIntervalVelocities() failed interval above botSurface " +
                                         botSurface.getName()  + " row " + row + " col " + col + ".", e, StsException.WARNING);
            return false;
        }
    }

    // not currently used:  when botDepthSurface is constructed, these intervalVelocity corrections are included
    // so they don't need to be recomputed/applied on db load.
/*
    private boolean correctIntervalVelocitiesToMarkers()
    {
        double botMarkerT, botMarkerZ;
         double topSurfaceT = 0, topSurfaceZ = 0;
         StsMarker marker = botTimeSurface.getMarker();
         if (marker == null)
         {
             StsException.systemError("StsVelocityGrid.correctIntervalVelocitiesToMarkers() failed." +
                                      " Surface " + botTimeSurface + " should have markers but doesn't.");
             return true;
         }

         Object[] wellMarkers = marker.getWellMarkers().getArrayList();
         int nWellMarkers = wellMarkers.length;
         if (wellMarkers == null || nWellMarkers == 0)
         {
             StsException.systemError("StsVelocityGrid.correctIntervalVelocitiesToMarkers() failed." +
                                      " Surface " + botTimeSurface + " should have markers but doesn't.");
             return true;
         }

         // For each marker, given its xy location, find the corresponding time
         // on the associated timeSurface.  This time corresponds to the depth
         // of the marker.  With the time and depth at the same xy of the surface
         // above, we can compute the interval velocity at that xy using the
         // the time and depth at the surface or datum just above.
         StsGridPoint gridPoint = new StsGridPoint(botTimeSurface);
         for (int n = 0; n < nWellMarkers; n++)
         {
             StsWellMarker wellMarker = (StsWellMarker) wellMarkers[n];
             if (wellMarker == null)
             {
                 continue;
             }
             StsPoint location = wellMarker.getLocation();
             float x = location.getX();
             float y = location.getY();
             float rowF = getRowCoor(y);
             float colF = getColCoor(x);
             float botSurfaceT = botTimeSurface.interpolateBilinearNoNulls(rowF, colF);
             botMarkerZ = location.getZ();
             topSurfaceT = getSurfaceTime(topTimeSurface, rowF, colF);
             topSurfaceZ = getSurfaceDepth(topDepthSurface, rowF, colF);
             float markerIntervalVelocity = (float) ( (botMarkerZ - topSurfaceZ) / (botSurfaceT - topSurfaceT));
             int row = (int) rowF;
             int col = (int) colF;
             double gridIntervalVelocity = StsMath.interpolateGridBilinear(velocities, rowF, colF);
             double adjustment = markerIntervalVelocity - gridIntervalVelocity;
             velocities[row][col] += adjustment;
             velocities[row][col + 1] += adjustment;
             velocities[row + 1][col + 1] += adjustment;
             velocities[row + 1][col] += adjustment;
             double adjustedGridVelocity = StsMath.interpolateGridBilinear(velocities, rowF, colF);
             System.out.println("DEBUG: Well marker " + marker.name + " well " + wellMarker.well.getName() +
                                " gridIntervalVelocity " + gridIntervalVelocity +
                                " adjustedGridIntervalVelocity " + adjustedGridVelocity);
        }
        return true;
    }
*/
    public boolean insideGrid(float rowF, float colF)
    {
        if (rowF < 0.0f || rowF > nRows - 1)
        {
            return false;
        }
        if (colF < 0.0f || colF > nCols - 1)
        {
            return false;
        }
        return true;
    }

    public float[] getVelocityRange()
    {
        return new float[]
            {
            velMin, velMax};
    }

    public double[][] getVelocities()
    {
        return velocities;
    }

    public double[] getRowGridlineValues(int nRow)
    {
        return velocities[nRow];
    }

    public double[] getColGridlineValues(int col)
    {
        double[] values = new double[nRows];
        for (int row = 0; row < nRows; row++)
        {
            values[row] = velocities[row][col];
        }
        return values;
    }
    /*
        public double[][] getErrorCorrections()
        {
            return errorCorrections;
        }
     */
}
