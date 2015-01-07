package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.UI.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

public class StsLogTrack extends StsObject implements ActionListener
{
	public StsObjectRefList logCurves;

	transient StsWellViewModel wellViewModel;
	//transient StsLogCurvesView logCurvesView = null;
	transient StsLogCurve selectedCurve;
	transient ButtonGroup selectCurveButtonGroup = new ButtonGroup();
	transient Color[] seismicColors = null;
	transient boolean displaySeismic = false;
	transient boolean displayTraces = false;
	transient boolean displaySynthetic = false;
	transient boolean displayDerived = false;
	transient String derivedLog = null;
	transient StsSeismicVolume seismicVolume = null;
	transient StsFloatList mdepthsList = null;
	transient StsObjectList colorsList = null;
	transient StsFloatList wiggleExcList = null;
	transient boolean colorscaleChanged = true;
	transient double windowMdepthMax, windowMdepthMin;
	transient StsLogCurve RHOB = null;
	transient StsLogCurve DT = null;
	transient int LTABLE = 8;
	transient int NTABLE = 513;
	static float table[][] = new float[513][8];
	static int tabled=0;
	transient String DTname=null;
	transient String RHOBname=null;
	/** maximum number of curves per track */
	static public int maxCurvesPerTrack = 5;

	public StsLogTrack()
	{
	}

	public StsLogTrack(StsWellViewModel wellViewModel, StsModel model, int nSubWindow)
	{
		super(false);

		this.wellViewModel = wellViewModel;
//		this.nTrack = nTrack;
		addToModel();
		logCurves = StsObjectRefList.constructor(maxCurvesPerTrack, 1, "logCurves", this);
		dbFieldChanged("logCurves", logCurves);
		StsLogCurve[] cv = getLogCurves();

		//logCurvesView = new StsLogCurvesView(wellViewModel, model, nSubWindow, this);
	}

	public boolean initialize(StsModel model)
	{
		return true;
	}

	public StsLogCurve[] getLogCurves()
	{
		if(logCurves.getSize() == 0)return new StsLogCurve[0];
		return(StsLogCurve[])logCurves.getCastList();
	}

	public StsLogCurve getSelectedLogCurve()
	{
		return selectedCurve;
	}

	public boolean removeLogCurve(StsLogCurve logCurve)
	{
		logCurves.delete(logCurve);
		if(logCurves.getSize() > 0)return false;
		currentModel.delete(this);
		return true;
	}

	public boolean removeAllLogCurves()
	{
		logCurves.deleteAll();
		currentModel.delete(this);
		return true;
	}

	public void setDisplayTraces(boolean b)
	{
		this.displayTraces = b;
	}

	public void setDisplaySynthetic(boolean b)
	{
		this.displaySynthetic = b;
	}

    public void setSyntheticNames(String DTname, String RHOBname)
	{
		this.DTname = DTname;
		this.RHOBname = RHOBname;
    }
	public void setDisplayDerived(boolean b, String derived)
	{
		this.displayDerived = b;
		this.derivedLog = derived;
	}

	/*
	 public void startGL()
	 {
	  logCurvesView.startGL();
	 }
	 */
	public void removeButtons()
	{
		selectCurveButtonGroup = null;
		selectCurveButtonGroup = new ButtonGroup();
	}

	public void addRadioButtonToGroup(JRadioButton radioButton)
	{
		selectCurveButtonGroup.add(radioButton);
		radioButton.addActionListener(this);

		/*
		 if(selectCurveButtonGroup.getButtonCount() == 1)
		 {
		  selectCurveButtonGroup.setSelected(radioButton.getModel(), true);
		  selectedCurve = logCurves[0];
		 }
		 */
	}

	public void initializeSelectedCurve()
	{
		StsLogCurve[] logCurveList = getLogCurves();
		if(logCurves.getSize() > 0)
		{

			if(logCurveList[0].getName().equals("Seismic"))
				displayTraces = true;
			if(logCurveList[0].getName().equals("Synthetic"))
				displaySynthetic = true;
		}
		if(logCurveList.length == 0)return;

		JRadioButton radioButton = (JRadioButton)selectCurveButtonGroup.getElements().nextElement();
		if(radioButton == null) return;

		selectCurveButtonGroup.setSelected(radioButton.getModel(), true);
		selectedCurve = (StsLogCurve)logCurves.getFirst();
	}

	public void actionPerformed(ActionEvent e)
	{
		Object source = e.getSource();
		if(source == null)return;

		if(source instanceof StsColorscale)
		{
			colorscaleChanged = true;
			wellViewModel.repaint();
		}
		else if(source instanceof JRadioButton)
		{
			if(logCurves.getSize() == 0)return;
			for(int n = 0; n < logCurves.getSize(); n++)
			{
				StsLogCurve logCurve = (StsLogCurve)logCurves.getElement(n);
				if(logCurve.getName() == e.getActionCommand())
				{
					selectedCurve = logCurve;
					wellViewModel.repaint();
					return;
				}
			}
		}
	}

	public double getValueFromPanelXFraction(double fraction)
	{
		if(selectedCurve == null)
		{
			return StsParameters.nullValue;
		}
		return selectedCurve.getValueFromPanelXFraction(fraction);
	}

	public boolean delete()
	{
		return super.delete();
	}

	public void seismicChanged()
	{
		colorscaleChanged = true;
	}

	public void drawSynthWiggle(GL gl, GLU glu, double[]mdepths, double[]vals, StsColor bg, StsColor fg)
	{

		// jbw display seismic


		StsWiggleDisplayProperties wiggleProperties = new StsWiggleDisplayProperties();
		wiggleProperties.setWiggleBackgroundColor(bg);
		wiggleProperties.setWigglePlusColor(fg);
        int nInterpolatedIntervals = 5;

        float[] interpolatedValues = StsTraceUtilities.computeCubicInterpolatedPoints(vals, nInterpolatedIntervals);
        float[] interpolatedMDepths = StsTraceUtilities.computeLinearInterpolatedPoints(mdepths, nInterpolatedIntervals);
		ortho2D(-1, 1, windowMdepthMax, windowMdepthMin, gl, glu);
		{
				//if(wiggleProperties.hasFill())
				StsTraceUtilities.drawFilledWiggleTraces(gl, interpolatedValues, interpolatedMDepths, 0.f, wiggleProperties); // horiz scale introduced above
				//if (wiggleProperties.getWiggleDrawLine())
				StsTraceUtilities.drawWiggleTraces(gl, interpolatedValues, interpolatedMDepths, 0.f, wiggleProperties); // horiz scale introduced above
			}
			gl.glDisable(GL.GL_LIGHTING);


		drawAxes(gl, -1, 1, true, 1.f, 1.f);
		popProjectionMatrix(gl);
	}

    public void drawWiggle(GL gl, GLU glu)
		{

			// jbw display seismic

			int nSeismicPoints = colorsList.getSize();
			float[] values = new float[nSeismicPoints];
            float[] zValues = new float[nSeismicPoints];
            //boolean wiggleSmoothCurve = true;
			for(int n = 0; n < nSeismicPoints; n++)
			{
				values[n] = wiggleExcList.getElement(n);
				zValues[n] = mdepthsList.getElement(n);
			}
			StsMath.normalizeAmplitude(values);
			StsWiggleDisplayProperties wiggleProperties = new StsWiggleDisplayProperties();
            StsProject project = currentModel.getProject();
            wiggleProperties.setWiggleBackgroundColor(project.getBackgroundColor());
			wiggleProperties.setWigglePlusColor(project.getForegroundColor());
			ortho2D(0., 7., windowMdepthMax, windowMdepthMin, gl, glu);
			{
				//if (wavePoints.length < 3) continue;

				//double[][] drawPoints;
				//if (wiggleSmoothCurve)
				//	drawPoints = StsTraceUtilities.computeInterpolatedDataPoints(wavePoints, pixelsPerInc, zInc, valueScale, valueOffset);
				//else
				//	drawPoints = StsTraceUtilities.getDisplayedData(wavePoints, displayedDataRange, valueScale, valueOffset);

				//                   float offset = (float)line.traceOffsets[nLineIndex];
				//if (adjustNMO) checkAdjustFlattening(drawPoints, nGatherIndex);

				//double[] traceMuteRange = getUnFlattenedMuteRange(nGatherIndex);
				//if (flatten) traceMuteRange = getFlattenedMuteRange(nGatherIndex);

				float x0;

				for(int iloop = 0; iloop < 6; iloop++)
				{
					x0 = iloop + 1;
					//if(wiggleProperties.hasFill())
					StsTraceUtilities.drawFilledWiggleTraces(gl, values, zValues, x0, wiggleProperties); // horiz scale introduced above
					//if (wiggleProperties.getWiggleDrawLine())
					StsTraceUtilities.drawWiggleTraces(gl, values, zValues, x0, wiggleProperties); // horiz scale introduced above
				}
				gl.glDisable(GL.GL_LIGHTING);
			}

			drawAxes(gl, 0.f, 7.f, true, 1.f, 1.f);
			popProjectionMatrix(gl);
		}

	public void drawLogTrack(GLAutoDrawable component, boolean displayValues, GL gl, GLU glu, StsMousePoint mousePoint,
                             int width, int height, StsWellViewModel wellViewModel, StsWellView wellView, ArrayList<StsPoint> drawPoints)
	{
		try
		{
			this.wellViewModel = wellViewModel;
			windowMdepthMax = wellViewModel.getWindowMdepthMax();
            windowMdepthMin = wellViewModel.getWindowMdepthMin(wellView);
			drawSeismic(component, gl, glu);
			drawSynthetic(component, gl, glu);
			drawDerivedCurves(component, gl, glu, displayValues, mousePoint, width, height);
			drawLogCurves(component, gl, glu, displayValues, mousePoint, width, height, drawPoints);
		}
		catch(Exception e)
		{
			StsException.outputException("drawLogTrack() failed.", e, StsException.WARNING);
		}
	}

	public static float geometricMean(float[] values, int start, int end)
			{
					int i;
					double product = 0.0;


					for(i=start; i<=end; i++)
							product *= values[i];

					return (float) Math.pow(product, (1/(end - start + 1)));
			}

    public static float[] applyGeoFilter(float[] vals, int n)
	{
		float[] ret = new float[vals.length];
		for (int i = n-1; i < vals.length; i++)
			ret[i]=geometricMean(vals,1,1+n-1);
		for (int i = 0; i < n-1; i++)
			ret[i]=vals[i];
		return vals;
    }

	private void drawSynthetic(GLAutoDrawable component, GL gl, GLU glu)
	{
		if (!displaySynthetic) return;

		if (DT == null || RHOB == null)
		{
		StsLogCurve[] logCurves = (StsLogCurve[])wellViewModel.well.getLogCurves().getCastListCopy();

		for (int i=0; i < logCurves.length; i++)
		{
			if (logCurves[i].getCurvename().equals(RHOBname))
				RHOB = logCurves[i];
			if (logCurves[i].getCurvename().equals(DTname))
			   DT = logCurves[i];
		}

		}
		if (RHOB == null)  return;
		float[] r   = RHOB.getValuesVectorFloats();
		if (r.length <= 0) return;
		float[] depths = RHOB.getDepthVectorFloats();
		float[] mdepths = RHOB.getMDepthVectorFloats();
		float [] dt = null;
		if (DT != null)
		   dt =  DT.getValuesVectorFloats();


		if (dt.length <= 0)
		{
			// use Gardner
			System.out.println("No DT curve, using basic Gardner");
			dt = new float[r.length];
			for (int i = 0; i < r.length; i++)
			{
				dt[i] = 1000000.f/(float)(Math.pow((r[i]/0.23),1./0.25));
			}
	   }
		if (r.length <= 0) return;
		float [] vels = new float[dt.length];
		for (int i = 0; i < vels.length; i++) // convert to ft/sec from microsec/ft
			vels[i] = 1000000.f /dt[i];

		// block logs. Should be Backus ?
		float b_vels [] = StsSeismicFilter.applyBoxFilter(vels, 1, vels.length, 8);
		float b_r [] =StsSeismicFilter.applyBoxFilter(r, 1, r.length, 8);

		// artifacts
		int last=vels.length-1;
		for (int i = 0; i < 4; i++)
		{
			b_vels[i] = vels[i];
			b_r[i] = r[i];
			b_vels[last-i] = vels[last-i];
			b_r[last-i] = r[last-i];
		}
		// min dt
		double tmax = 0;
		double dtmin = 1.E30;
		for (int i = 1; i < depths.length; i++)
		{
			double tempdt = 0.0;
			{
			if (b_vels[i-1] < 1.e30) // missing
				tempdt =  2.0 * (depths[i]-depths[i-1])/b_vels[i-1];
				dtmin = dtmin < tempdt ? dtmin: tempdt;
				tmax += tempdt;
			}

		}
		if (dtmin < 0.0001) dtmin = 0.0001; // effectively limits hi freq
		int ntmax = (int)( 2.5 + tmax/dtmin);
		double ttotal = 0.;
		//double maxtrace[] = new double[ntmax];
		double rc[] = new double[depths.length];
		float  rctimes[] = new float[depths.length];

		float dtout = 0.04f;
		int ntout = (int) (0.5f + tmax/dtout);


		for (int i = 1; i < depths.length; i++)
		{
			if (b_vels[i-1] < 1.e30)
			{
			/*
				ttotal += 2.0 * (depths[i] - depths[i - 1]) / b_vels[i - 1];
				int itempt = (int)((double)ntmax * ttotal / (double)tmax);
				if(itempt >= ntmax)
				{
					itempt = ntmax - 1;
					//System.out.println("over " + ttotal);
				}
				else
				{
					//System.out.println("i "+i+" ttotal "+ttotal);
					maxtrace[itempt] += (b_vels[i] * b_r[i] - b_vels[i - 1] * b_r[i - 1]) / (b_vels[i] * b_r[i] + b_vels[i - 1] * b_r[i - 1]);
				}
			*/
		         rc[i] = (b_vels[i] * b_r[i] - b_vels[i - 1] * b_r[i - 1]) / (b_vels[i] * b_r[i] + b_vels[i - 1] * b_r[i - 1]);
				 rctimes[i] = wellViewModel.well.getTimeFromMDepth(mdepths[i]);
			}
		}

		// rctimes are irregular map back to standard sampling
		float t0 = wellViewModel.well.getTimeFromMDepth(mdepths[0]);
		rctimes[0]=t0;
		float eventimes[] = new float[mdepths.length];
		double maxtrace[] = new double[mdepths.length];
		for (int i = 0; i < mdepths.length; i++)
			 maxtrace[i]=0;
		float evendt = (rctimes[mdepths.length-1]-rctimes[0])/(float)(mdepths.length -1);
		for (int i = 0; i < mdepths.length; i++)
		{
			eventimes[i] = i * evendt;
			int timeind = (int)((rctimes[i] - t0)/evendt + .5f);
			if (timeind >= mdepths.length) timeind = mdepths.length-1;
			if (timeind < 0) timeind = 0;
			maxtrace[timeind] += rc[i];
		}
		//StsSeismicFilter.applyBoxFilter(maxtrace, 1, maxtrace.length, 8);
		double [] ricker = StsSeismicFilter.ricker(60.,(float)(evendt/1000.f));
		double [] newdepths = new double[maxtrace.length];
		for (int i = 0; i < maxtrace.length; i++)
			newdepths[i] = mdepths[i];


		double [] filtered= new double[maxtrace.length];
		StsConvolve.convolve(maxtrace, ricker, filtered, maxtrace.length, 0, maxtrace.length-1, ricker.length, ricker.length/2);
        StsProject project = currentModel.getProject();
		drawSynthWiggle( gl,  glu, newdepths, filtered, project.getBackgroundColor(), project.getForegroundColor());
		//drawSynthWiggle( gl,  glu, newdepths, maxtrace, Color.WHITE, Color.RED);
	}

	private void drawSynthetic2(GLAutoDrawable component, GL gl, GLU glu)
	{
		if (!displaySynthetic) return;

	    if (DT == null || RHOB == null)
		{
		StsLogCurve[] logCurves = (StsLogCurve[])wellViewModel.well.getLogCurves().getCastListCopy();

		for (int i=0; i < logCurves.length; i++)
		{
			if (logCurves[i].getCurvename().equals(RHOBname))
			    RHOB = logCurves[i];
			if (logCurves[i].getCurvename().equals(DTname))
			   DT = logCurves[i];
		}

		}
		if (RHOB == null)  return;
		float[] r   = RHOB.getValuesVectorFloats();
		if (r.length <= 0) return;
		float[] depths = RHOB.getDepthVectorFloats();
		float[] mdepths = RHOB.getMDepthVectorFloats();
		float [] dt = null;
		if (DT != null)
		   dt =  DT.getValuesVectorFloats();


		if (dt.length <= 0)
		{
			// use Gardner
			System.out.println("No DT curve, using basic Gardner");
			dt = new float[r.length];
			for (int i = 0; i < r.length; i++)
			{
				dt[i] = 1000000.f/(float)(Math.pow((r[i]/0.23),1./0.25));
			}
	   }
		if (r.length <= 0) return;
		float [] vels = new float[dt.length];
		for (int i = 0; i < vels.length; i++) // convert to ft/sec from microsec/ft
			vels[i] = 1000000.f /dt[i];

	    // block logs. Should be Backus ?
		float b_vels [] = StsSeismicFilter.applyBoxFilter(vels, 1, vels.length, 8);
		float b_r [] =StsSeismicFilter.applyBoxFilter(r, 1, r.length, 8);
		//float b_vels [] = applyGeoFilter(vels, 8);
		//float b_r [] =applyGeoFilter(r, 8);

        // artifacts
		int last=vels.length-1;
		for (int i = 0; i < 4; i++)
		{
			b_vels[i] = vels[i];
			b_r[i] = r[i];
			b_vels[last-i] = vels[last-i];
			b_r[last-i] = r[last-i];
		}
		// min dt
		double tmax = 0;
		double dtmin = 1.E30;

			for(int i = 1; i < depths.length; i++)
			{

				double tempdt = 0.0;
				{
					if(b_vels[i - 1] < 1.e30) // missing
						tempdt = 2.0 * (depths[i] - depths[i - 1]) / b_vels[i];
					dtmin = dtmin < tempdt ? dtmin : tempdt;
                    tmax += tempdt;
				}

			}

			if(dtmin < 0.0001)dtmin = 0.0001; // effectively limits hi freq
			int ntmax = (int)(2.5 + tmax / dtmin);
			double ttotal = 0.;
			double maxtrace[] = new double[ntmax];
			for(int i = 0; i < ntmax; i++)
				maxtrace[i] = 0.f;



/* su style */



		//float dtout = 0.01f;
		//int ntout = (int) (0.5f + tmax/dtout);

		float t0 = wellViewModel.well.getTimeFromMDepth(mdepths[0]);
		ttotal = 0;
		for (int i = 1; i < depths.length; i++)
		{
			if (b_vels[i-1] < 1.e30)
			{
				ttotal += 2.0 * (depths[i] - depths[i - 1]) / (b_vels[i]);
				int itempt = (int)((double)ntmax * ttotal / (double)tmax);
				if(itempt >= ntmax)
				 itempt = ntmax - 1;
				else
				 maxtrace[itempt] += (b_vels[i] * b_r[i] - b_vels[i - 1] * b_r[i - 1]) / (b_vels[i] * b_r[i] + b_vels[i - 1] * b_r[i - 1]);
			}
		}



/*
	 float tmp1 = b_r[0] * b_vels[0];
	 float[] t = new float[b_r.length];
	 double[] rc = new double[b_r.length];
	 rc[0]=0;
	 t[0]=0;
	 float t0 = wellViewModel.well.getTimeFromMDepth(mdepths[0]);
	 for (int i = 1; i < depths.length; i++)
	  {
		  if (b_vels[i-1] < 1.e30)
		  {
			  float tmp2 = b_r[i] * b_vels[i];
			  rc[i] = (tmp2-tmp1)/(tmp2+tmp1);
			  tmp1=tmp2;
			  t[i]=t0+t[i-1]+2.f*(depths[i] - depths[i - 1]) / b_vels[i];
			  System.out.println("rc= "+rc[i]+" t[i]= "+1000.f*t[i]);
	      }

	  }

	  float tn = wellViewModel.well.getTimeFromMDepth(mdepths[mdepths.length-1]);
	  float dt2=0.002f;
	  int nt = (int) ((tn - t0)/dt2 + .5);
	  float[] tout = new float[nt];
	  float[] aux_trace = new float[nt];
	  // create array of equisampled times for output trace
	  for (int k = 0; k < nt; k++)
		  tout[k]=t0+(k*dt2);


	  //interpolate times to uniform sampling interval
        //com.Sts.Utilities.SincInterpolate.ints8r(depths.length,dt2,0.f,t,rc,0.0f,0.0f,nt,tout,aux_trace);
*/


        double [] ricker = StsSeismicFilter.ricker(60.,(float)dtmin);
        double [] newdepths = new double[maxtrace.length];
		double [] filtered= new double[maxtrace.length];
		StsConvolve.convolve(maxtrace, ricker, filtered, maxtrace.length, 0, maxtrace.length-1, ricker.length, ricker.length/2);

		for (int i =0; i < maxtrace.length; i++)
		{
			newdepths[i] = wellViewModel.well.getMDepthFromTime((float) (t0+ (1000.f*i*dtmin)));
		}


        StsMath.normalizeAmplitude(filtered);
        StsProject project = currentModel.getProject();
        drawSynthWiggle( gl,  glu, newdepths, filtered, project.getBackgroundColor(), project.getForegroundColor());
		//drawSynthWiggle( gl,  glu, newdepths, maxtrace, Color.WHITE, Color.RED);
	}
	private void drawSynthetic3(GLAutoDrawable component, GL gl, GLU glu)
		{
			if (!displaySynthetic) return;


			if (DT == null || RHOB == null)
			{
			StsLogCurve[] logCurves = (StsLogCurve[])wellViewModel.well.getLogCurves().getCastListCopy();

			for (int i=0; i < logCurves.length; i++)
			{
				if (logCurves[i].getCurvename().equals(RHOBname))
					RHOB = logCurves[i];
				if (logCurves[i].getCurvename().equals(DTname))
				   DT = logCurves[i];
			}

			}
			if (DT == null || RHOB == null)  return;
			float[] depths = RHOB.getDepthVectorFloats();
			float[] mdepths = RHOB.getMDepthVectorFloats();
			float[] vel =  DT.getValuesVectorFloats();
			float[] rho   = RHOB.getValuesVectorFloats();
			float[] deltaz = new float[depths.length];
			// min dt
			double tmax = 0;
			double dtmin = 1.E30;
			for (int i = 1; i < depths.length; i++)
			{
				deltaz[i]=depths[i]-depths[i-1];
			}
			// compute # samp in seismogram
			float dt = 0.004f;
			tmax = 2.0*deltaz[1]/vel[1];
			for (int i = 1; i < depths.length; i++)
				tmax=tmax + 2.0 * deltaz[i]/vel[i];

		    int nt =  (int) (tmax/dt);
			System.out.println("nt is "+nt);
			float[] trace=new float[nt];
			int nlyr=nt;
			float freq=25.f;
			createSynthetic(nt,nlyr,freq, dt,vel,rho,deltaz,trace);

			//drawSynthWiggle( gl,  glu, mdepths, trace);
	}
	private void createSynthetic(int nt, int nlyr, float freq, float dt, float[] vel, float[] rho,
								 float[]deltaz, float[]trace)
	{
		float[] rc = new float[nlyr];
		float [] t = new float[nlyr];
		float[] tout = new float[nt];
		float [] aux_trace = new float[nt];
		float[] wavelet = null;

		float tmp1, tmp2;
		tmp1=rho[0]*vel[0];
		rc[0]=0.0f;
		t[0]=2.f*deltaz[0]/vel[0];

		for (int i = 1; i < nlyr; i++)
		{
			tmp2 = rho[i]*vel[i];
			rc[i]=(tmp2-tmp1)/(tmp2+tmp1);
			tmp1=tmp2;
			t[i]=t[i-1]+2.0f*deltaz[i]/vel[i];
		}

		// 8 pt sinc interp
		//interp8SU(nlyr,dt,t,rc,0.0f, 0.0f,nt,tout,aux_trace);

		// ricker
		int lw, iw;
		//wavelet = ruckerSU(freq,dt,lw,iw);

		// convolve
		//trace= convolveSU(lw,0,wavelet,nt,0,aux_trace,nt,0);

	}



	private double[] convolveSU(int lx, int ifx, float[]x, int ly, int ify, float[]y, int lz, int ifz)
	{
		double[]z = new double[lz];
		int ilx = ifx + lx -1;
		int ily = ify + ly -1;
		int ilz = ifz + lz -1;
		float zsum= 0.f;

		for (int i = 0; i < lz; i++)
		{
			int jlow = i - ily;
			if (jlow < ifx) jlow=ifx;
			int jhigh = i-ify;
			if (jhigh > ilx) jhigh = ilx;
			zsum = 0.0f;
			for (int j = jlow; j <=jhigh; ++j)
				zsum += x[j]*y[i-j];

		     z[i]=zsum;
		}

		return z;
	}
	private void drawSeismic(GLAutoDrawable component, GL gl, GLU glu)
	{

		if(!displaySeismic && !displayTraces)return;
		StsSeismicVolume currentSeismicVolume = (StsSeismicVolume)currentModel.getCurrentObject(StsSeismicVolume.class);
		if(currentSeismicVolume == null)return;

		StsPoint[] points = wellViewModel.well.getRotatedPoints();
		int nPoints = points.length;

		if(seismicVolume != currentSeismicVolume)
		{
			//if (line2d != null) line2d.removeActionListener(this);
			if(currentSeismicVolume != null)currentSeismicVolume.getColorscale().addActionListener(this);
			seismicVolume = currentSeismicVolume;
			colorscaleChanged = true;
		}
		//if(colorscaleChanged)
		// jbw until we work out callbacks for editing, do this every time
		if (true)
		{
			float tMin = seismicVolume.getZMin();
			float tMax = seismicVolume.getZMax();
			float tInc = seismicVolume.getZInc();
			float t1 = points[0].getT();
			float m1 = points[0].getM();
			float x1 = points[0].getX();
			float y1 = points[0].getY();
			float f1 = (t1 - tMin) / tInc;
			boolean out1 = t1 < tMin || t1 > tMax;
			int i0, i1;
			mdepthsList = new StsFloatList(100, 100);
			colorsList = new StsObjectList(100, 100);
			wiggleExcList = new StsFloatList(100, 100);
			for(int n = 1; n < nPoints; n++)
			{
				float t0 = t1;
				float f0 = f1;
				float m0 = m1;
				float x0 = x1;
				float y0 = y1;
				boolean out0 = out1;
				t1 = points[n].getT();
				m1 = points[n].getM();
				x1 = points[n].getX();
				y1 = points[n].getY();
				f1 = (t1 - tMin) / tInc;

				out1 = t1 < tMin || t1 > tMax;
				if(out0 && out1)continue;

				if(f1 > f0)
				{
					i0 = StsMath.ceiling(f0);
					i1 = StsMath.below(f1);
					if(i1 < i0)continue;

					float f = (i0 - f0) / (f1 - f0);
					float m = f * (m1 - m0) + m0;
					float x = f * (x1 - x0) + x0;
					float y = f * (y1 - y0) + y0;
					float t = i0 * tInc + tMin;
					float df = 1 / (f1 - f0);
					float dm = df * (m1 - m0);
					float dx = df * (x1 - x0);
					float dy = df * (y1 - y0);
					for(int i = i0; i <= i1; i++, t += tInc, m += dm, x += dx, y += dy)
					{
						Color color = seismicVolume.getColor(x, y, t);
						colorsList.add(color);
						mdepthsList.add(m);
						wiggleExcList.add(3.0f * (float)seismicVolume.getValue(x, y, t) / 255.f);
					}
				}
				else if(f1 < f0)
				{
					i0 = StsMath.floor(f0);
					i1 = StsMath.above(f1);
					if(i1 > i0)continue;

					float f = (i0 - f0) / (f1 - f0);
					float m = f * (m1 - m0) + m0;
					float x = f * (x1 - x0) + x0;
					float y = f * (y1 - y0) + y0;
					float t = i0 * tInc;
					float df = -1 / (f1 - f0);
					float dm = df * (m1 - m0);
					float dx = df * (x1 - x0);
					float dy = df * (y1 - y0);

					for(int i = i0; i >= i1; i--, t -= tInc, m += dm, x += dx, y += dy)
					{
						Color color = seismicVolume.getColor(x, y, t);
						colorsList.add(color);
						mdepthsList.add(m);
						wiggleExcList.add(3.0f * (float)seismicVolume.getValue(x, y, t) / 255.f);
					}
				}
			}
			colorscaleChanged = false;
		}

		if(displaySeismic)
		{
			ortho2D(0.0, 1.0, windowMdepthMax, windowMdepthMin, gl, glu);

			float m1 = mdepthsList.getElement(0);
			float m = m1;
			gl.glBegin(GL.GL_QUAD_STRIP);
			int nSeismicPoints = colorsList.getSize();
			if(true)
				for(int n = 0; n < nSeismicPoints; n++)
				{
					Color color = (Color)colorsList.getElement(n);
					StsColor.setGLJavaColor(gl, color);
					gl.glVertex2f(0.0f, m);
					gl.glVertex2f(1.0f, m);
					if(n < nSeismicPoints - 1)
					{
						float m0 = m1;
						m1 = mdepthsList.getElement(n + 1);
                        m = (m0 + m1) / 2;
					}
					else
						m = m1;
					gl.glVertex2f(0.0f, m);
					gl.glVertex2f(1.0f, m);
				}
			gl.glEnd();
			popProjectionMatrix(gl);
		}
		drawWiggle(gl, glu);
	}

	private void drawAxes(GL gl,float gridMin, float gridMax, boolean isLinear, float boldInc, float gridInc)
	{
	double unitsPerLabel = wellViewModel.zoomLevel.unitsPerLabel;
	double firstLabeledZ = StsMath.intervalRoundUp(windowMdepthMin, unitsPerLabel);
	double lastLabeledZ = StsMath.intervalRoundDown(windowMdepthMax, unitsPerLabel);
	int nGridValues = (int) (1 + (lastLabeledZ - firstLabeledZ) / unitsPerLabel);
	double value;
	// draw horizontal gridLines across at each labeled depth value
	// StsColor c = StsColor.getInverseStsColor(currentModel.project.getForegroundColor());
	//StsColor.BLACK.setGLColor(gl);
    currentModel.getProject().getForegroundColor().setGLColor(gl);
    //c.setGLColor(gl);
	gl.glBegin(GL.GL_LINES);
	value = firstLabeledZ;
	if (isLinear)
	{
		for (int n = 0; n < nGridValues; n++)
		{
			gl.glVertex2d(gridMin, value);
			gl.glVertex2d(gridMax, value);
			value += unitsPerLabel;
		}
	}
	/*
	else
	{
		for (int n = 0; n < nGridValues; n++)
		{
			gl.glVertex2d(logGridMin, value);
			gl.glVertex2d(logGridMax, value);
			value += unitsPerLabel;
		}
	}
	*/
	gl.glEnd();
	// for currently selectedCurve, draw vertical gridLines down
	if (true)
	{
		if (isLinear)
		{
			value = gridMin;
			nGridValues = (int) (1 + (gridMax - gridMin) / gridInc);
			double firstBoldValue = StsMath.intervalRoundUp(gridMin, boldInc);
			int nextBoldLine = (int) ( (firstBoldValue - gridMin) / gridInc);
			boldInc = (int) (boldInc / gridInc);
			gl.glBegin(GL.GL_LINES);
			for (int n = 0; n < nGridValues; n++)
			{
				if (n == nextBoldLine)
				{

                    //StsColor c2 = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
					//c2.setGLColor(gl);
                    currentModel.getProject().getForegroundColor().setGLColor(gl);
                    gl.glLineWidth(3.f);
					nextBoldLine += boldInc;
				}
				else
				{
					gl.glLineWidth(1.f);
					StsColor.GREY.setGLColor(gl);
				}
				gl.glVertex2d(value, windowMdepthMin);
				gl.glVertex2d(value, windowMdepthMax);
				value += gridInc;
			}
			gl.glEnd();
		}
		/*
		else
		{
			nGridValues = (int) (1 + (logGridMax - logGridMin));
			value = logGridMin;
			double drawValue;
			gl.glBegin(GL.GL_LINES);
			for (int n = 0; n < nGridValues; n++)
			{
				StsColor.WHITE.setGLColor(gl);
				gl.glVertex2d(value, windowMdepthMin);
				gl.glVertex2d(value, windowMdepthMax);
				StsColor.GREY.setGLColor(gl);
				drawValue = value + StsMath.log10(5);
				gl.glVertex2d(drawValue, windowMdepthMin);
				gl.glVertex2d(drawValue, windowMdepthMax);
				value += 1;
			}
			gl.glEnd();
		}
		*/
	}
	}
	private void drawDerivedCurves(GLAutoDrawable component, GL gl, GLU glu, boolean displayValues, StsMousePoint mousePoint, int width, int height)
	{
		if (!displayDerived) return;
		float gridMin = 0.f;
		float gridMax = 1.f;
		float boldInc = 1.f;
		float gridInc = 1.f;
		// jbw recompute rotated points
		//StsPoint[] points = wellViewModel.well.getPointsFromLineVertices();
		StsPoint[] points = wellViewModel.well.getRotatedPoints();
		int nPoints = points.length;
		float t0 = points[0].getT();
        float m0 = points[0].getM();
		float t1 = points[nPoints-1].getT();
		if (t0 < 0)
			t0 = 0;
		t1 = t1 * 1.1f;
		// assume monotonic;


		if (derivedLog.equals("TIME/DEPTH"))
		{
			ortho2D (t0, t1, windowMdepthMax, windowMdepthMin, gl, glu);
			gridMin = t0;
			gridMax = t1;
			gridInc = (t1-t0)/10.f;
			boldInc = gridInc/3.f;

			//StsColor c = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
			//c.setGLColor(gl);
            currentModel.getProject().getForegroundColor().setGLColor(gl);
            //StsColor.setGLJavaColor(gl, Color.BLACK);
			gl.glBegin(GL.GL_LINE_STRIP);
			for(int n = 0; n < nPoints; n++)
			{
				gl.glVertex2f(points[n].getT(), points[n].getM());
			}
			gl.glEnd();
			drawAxes(gl, gridMin, gridMax, true, boldInc, gridInc);
		}

		if (derivedLog.equals("Interval Velocity"))
		{
			ortho2D(0.f, 10000.f, windowMdepthMax, windowMdepthMin, gl, glu);
			gridMin=0.f;
			gridMax=20000.f;
			boldInc = 2000.f;
			gridInc = 1000.f;

	/*
			StsLogCurve[] logCurves = (StsLogCurve[])wellViewModel.well.getLogCurves().getCastListCopy();
			StsLogCurve TIME = null;
			for (int i=0; i < logCurves.length; i++)
			{
				if (logCurves[i].getCurvename().equals("TIME"))
					TIME = logCurves[i];
			}

	        if (TIME == null)  return;
            float[] r   = TIME.getValuesVectorFloats();
            if (r.length <= 0) return;
			float[] d = TIME.getDepthFloatVector().getValues();

			StsColor.setGLJavaColor(gl, Color.BLACK);
			gl.glBegin(GL.GL_LINE_STRIP);
			for (int n = 1; n < r.length; n++)
			{
				//float dt = points[n].getT() - points[n-1].getT();
				float dt = r[n] - r[n-1];
				float dz = d[n] - d[n-1];
				float iv;
				if (dt > 0)
				   iv = 1000.f*dz/dt;
			    else
				   iv = 20000.f;
				gl.glVertex2f(iv,points[n].getM());
			}
	*/
            //StsColor c = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
            //c.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 1; n < nPoints; n++)
            {
	           float dt = points[n].getT() - points[n-1].getT();
	           float dz = points[n].getZ() - points[n-1].getZ();
	           float iv = 1000.f*dz/dt;
	           gl.glVertex2f(iv,points[n].getM());
            }
            gl.glEnd();

			drawAxes(gl,gridMin, gridMax, true, boldInc, gridInc);


		}

		popProjectionMatrix(gl);

		if (displayValues)
		{
			wellViewModel.push2DOrtho(gl, glu, width, height);
			displayValues(gl, gridMin, gridMax, true, mousePoint, width, height);
		}
		else
		{
			wellViewModel.push2DOrtho(gl, glu, width, height);
			displayValues2(gl, gridMin, gridMax, true, mousePoint, width, height);
		}
		wellViewModel.pop2DOrtho(gl);
	}

	private void drawLogCurves(GLAutoDrawable component, GL gl, GLU glu, boolean displayValues, StsMousePoint mousePoint, int width, int height, ArrayList<StsPoint> drawPoints)
	{
		// draw each of the log curves; for selected
		gl.glLineWidth(1.0f);
		StsLogCurve[] logCurves = getLogCurves();
		StsLogCurve selectedLogCurve = getSelectedLogCurve();
		for (int nCurve = 0; nCurve < logCurves.length; nCurve++)
		{
			StsLogCurve logCurve = logCurves[nCurve];
			if (logCurve == null)
			{
				continue;
			}
			float[] mdepths = getMDepthVectorFloats(logCurve);
			if (mdepths == null)
			{
				float[] depths = logCurve.getDepthVectorFloats();
				if(depths == null) continue;
				mdepths = wellViewModel.well.getMDepthsFromDepths(depths);
				if(mdepths == null) continue;
			}
            float[] values = logCurve.getValuesVectorFloats();
			if (values == null)continue;

			StsLogCurveType curveType = logCurve.getLogCurveType();
			double[] scale = curveType.getScale();
			double gridMin = scale[0];
			double gridMax = scale[1];
			double gridInc = scale[2];
			double boldInc = scale[3];
			boolean isLinear = curveType.isLinear();
			double logGridMin = 0, logGridMax = 1;

			if (isLinear)
			{
				ortho2D(gridMin, gridMax, windowMdepthMax, windowMdepthMin, gl, glu);
			}
			else
			{
				logGridMin = StsMath.log10(gridMin);
				logGridMax = StsMath.log10(gridMax);
				ortho2D(logGridMin, logGridMax, windowMdepthMax, windowMdepthMin, gl, glu);
			}
			int nValues = values.length;
			logCurve.getStsColor().setGLColor(gl);
			boolean isDrawing = false;
			for (int n = 0; n < nValues; n++)
			{
				if (values[n] == StsParameters.nullValue)
				{
					if (isDrawing)
                    {
                        isDrawing = false;
                        gl.glEnd();
                    }
					continue;
				}
				else if (!isDrawing)
				{
					isDrawing = true;
					gl.glBegin(GL.GL_LINE_STRIP);
				}
				if (isLinear)
				{
					gl.glVertex2f(values[n], mdepths[n]);
				}
				else
				{
					gl.glVertex2d(StsMath.log10(values[n]), (double) mdepths[n]);
				}
			}
			if (isDrawing)
			{
				gl.glEnd();
			}
			double unitsPerLabel = wellViewModel.zoomLevel.unitsPerLabel;
			double firstLabeledZ = StsMath.intervalRoundUp(windowMdepthMin, unitsPerLabel);
			double lastLabeledZ = StsMath.intervalRoundDown(windowMdepthMax, unitsPerLabel);
			int nGridValues = (int) (1 + (lastLabeledZ - firstLabeledZ) / unitsPerLabel);
			double value;
			// draw horizontal gridLines across at each labeled depth value
			StsColor.GREY.setGLColor(gl);
			gl.glBegin(GL.GL_LINES);
			value = firstLabeledZ;
			if (isLinear)
			{
				for (int n = 0; n < nGridValues; n++)
				{
					gl.glVertex2d(gridMin, value);
					gl.glVertex2d(gridMax, value);
					value += unitsPerLabel;
				}
			}
			else
			{
				for (int n = 0; n < nGridValues; n++)
				{
					gl.glVertex2d(logGridMin, value);
					gl.glVertex2d(logGridMax, value);
					value += unitsPerLabel;
				}
			}
			gl.glEnd();
			// for currently selectedCurve, draw vertical gridLines down
			if (logCurve == getSelectedLogCurve())
			{
				if (isLinear)
				{
					value = gridMin;
					nGridValues = (int) (1 + (gridMax - gridMin) / gridInc);
					double firstBoldValue = StsMath.intervalRoundUp(gridMin, boldInc);
					int nextBoldLine = (int) ( (firstBoldValue - gridMin) / gridInc);
					boldInc = (int) (boldInc / gridInc);
					gl.glBegin(GL.GL_LINES);
                    StsColor.GREY.setGLColor(gl);
					for (int n = 0; n < nGridValues; n++)
					{
						if (n == nextBoldLine)
						{
						//	StsColor c = StsColor.getInverseStsColor(currentModel.getProject().getForegroundColor());

						//	c.setGLColor(gl);
							nextBoldLine += boldInc;
						}
                        /*
						else
						{
							StsColor.GREY.setGLColor(gl);
						}
						*/
						gl.glVertex2d(value, windowMdepthMin);
						gl.glVertex2d(value, windowMdepthMax);
						value += gridInc;
					}
					gl.glEnd();
				}
				else
				{
					nGridValues = (int) (1 + (logGridMax - logGridMin));
					value = logGridMin;
					double drawValue;
					gl.glBegin(GL.GL_LINES);
					for (int n = 0; n < nGridValues; n++)
					{
						//StsColor c = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
						//c.setGLColor(gl);
                        currentModel.getProject().getForegroundColor().setGLColor(gl);
                        gl.glVertex2d(value, windowMdepthMin);
						gl.glVertex2d(value, windowMdepthMax);
						StsColor.GREY.setGLColor(gl);
						drawValue = value + StsMath.log10(5);
						gl.glVertex2d(drawValue, windowMdepthMin);
						gl.glVertex2d(drawValue, windowMdepthMax);
						value += 1;
					}
					gl.glEnd();
				}

                if(drawPoints != null && drawPoints.size() > 0)
                {
                    StsColor.WHITE.setGLColor(gl);
                    gl.glPointSize(4.0f);
                    gl.glBegin(GL.GL_POINTS);
                    for(StsPoint drawPoint : drawPoints)
                        gl.glVertex2fv(drawPoint.v, 0);
                    gl.glEnd();
                }

				popProjectionMatrix(gl);

				if (displayValues)
				{
					wellViewModel.push2DOrtho(gl, glu, width, height);
					displayValues(gl, gridMin, gridMax, isLinear, mousePoint, width, height);
				}
				else
				{
					wellViewModel.push2DOrtho(gl, glu, width, height);
					displayValues2(gl, gridMin, gridMax, isLinear, mousePoint, width, height);
				}

				wellViewModel.pop2DOrtho(gl);
			}
			else
				popProjectionMatrix(gl);
		}

	}

	private void ortho2D(double left, double right, double bottom, double top, GL gl, GLU glu)
	{
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPushMatrix();
		gl.glLoadIdentity();
		glu.gluOrtho2D(left, right, bottom, top);
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

	private void popProjectionMatrix(GL gl)
	{
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPopMatrix();
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

	public float[] getMDepthVectorFloats(StsLogCurve curve)
	{
		if (curve.getMDepthVector() == null)
			return null;

		// Already Loaded and converted to project units
		StsFloatVector floatVector = curve.getMDepthVector().getValues();
		if(floatVector != null && floatVector.getValues() != null)
			return floatVector.getValues();

		// Load and convert to project units
//          String binaryDataDir = currentModel.getProject().getBinaryFullDirString();
		String binaryDataDir = curve.getStsDirectory() + currentModel.getProject().getBinaryDirString();
		curve.getMDepthVector().readBinaryFile(binaryDataDir, true);

		float scalar = currentModel.getProject().getDepthScalar(curve.getMDepthVector().getUnits());
		curve.getMDepthVector().convertToProjectUnits(scalar);

		floatVector = curve.getMDepthVector().getValues();
		if(floatVector == null)
			return null;
		return floatVector.getValues();
	}

	private void displayValues(GL gl, double gridMin, double gridMax, boolean isLinear, StsMousePoint mousePoint, int width, int height)
	{
		if (!wellViewModel.cursorPicked) return;
		//StsLogCurvesGLPanel displayPanel = getLogTrackView().getWellCurveDisplayPanel();
		//StsMousePoint mousePoint = displayPanel.getMousePoint();

		int glMouseX = mousePoint.x;
		int glMouseY = wellViewModel.cursorY;
		int locationX = StsMath.minMax(glMouseX, 0, Math.max(0, width - 70));
		int locationY = StsMath.minMax(glMouseY, 2, Math.max(0, height - 10));
		double depth = wellViewModel.getMdepthFromGLMouseY(glMouseY);
		double fraction = (double) glMouseX / width;
		double value = getValueFromPanelXFraction(fraction);
		if (value == StsParameters.nullValue)
			value = (fraction * (gridMax-gridMin)) + gridMin;
		//String depthLabel = wellViewModel.labelFormatter.format(depth);
		String depthLabel = wellViewModel.getZStringFromGLMouseY(glMouseY);
		String valueLabel = wellViewModel.labelFormatter.format(value);
		String mouseXLabel = wellViewModel.labelFormatter.format(glMouseX);
		String mouseYLabel = wellViewModel.labelFormatter.format(glMouseY);
        currentModel.getProject().getForegroundColor().setGLColor(gl);
        //StsColor c = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
		//c.setGLColor(gl);
		int mouseGLy = locationY;
		StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy, depthLabel + "   " + valueLabel);
//		StsGLDraw.highlightedFontHelvetica12(gl, locationX, mouseGLy, depthLabel + " " + valueLabel);
		//		StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy, mouseXLabel + " " + mouseYLabel);
		//		StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy-20, depthLabel + " " + valueLabel);
	}

    private void displayValues2(GL gl, double gridMin, double gridMax, boolean isLinear, StsMousePoint mousePoint, int width, int height)
	{
		if (!wellViewModel.cursorPicked) return;

		int mouseX = width/3;
		int mouseY = wellViewModel.cursorY;
		int locationX = StsMath.minMax(mouseX, 0, Math.max(0, width - 100));
		int locationY = StsMath.minMax(mouseY, 2, Math.max(0, height - 10));
		//double depth = wellViewModel.getZFromCursor();
		double fraction = (double) mouseX / width;
		double value = getValueFromPanelXFraction(fraction);
		//String depthLabel = wellViewModel.labelFormatter.format(depth);
		String depthLabel = wellViewModel.getMdepthStringFromGLCursorY();
		String valueLabel = wellViewModel.labelFormatter.format(value);
		// String mouseXLabel = wellViewModel.labelFormatter.format(mouseX);
		// String mouseYLabel = wellViewModel.labelFormatter.format(glMouseY);
        currentModel.getProject().getForegroundColor().setGLColor(gl);
        // StsColor c = StsColor.getInverseStsColor(currentModel.project.getStsWellPanelColor());
		// c.setGLColor(gl);
		// int mouseGLy = locationY;

//		StsGLDraw.fontHelvetica12(gl, value, depth, depthLabel + " " + valueLabel);
		StsGLDraw.fontHelvetica12(gl, locationX, locationY, depthLabel + "   " + valueLabel);
		//		StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy, mouseXLabel + " " + mouseYLabel);
		//		StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy-20, depthLabel + " " + valueLabel);
	}

    public boolean isDisplayingLogCurve(StsLogCurve logCurve)
    {
        return logCurves.hasObject(logCurve);
    }

	public void itemStateChanged(ItemEvent e)
	{
		if (e.getItem() instanceof StsColorscale)
		{
			colorscaleChanged = true;
			wellViewModel.repaint();
		}
	}
}
