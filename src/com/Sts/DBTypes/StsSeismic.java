package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.IO.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.nio.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public abstract class StsSeismic extends StsSeismicBoundingBox
{
	/** All the colorscales being used */
	protected StsObjectRefList colorscales = null;
	/** The current colorscale viewed on the object panel */
	protected StsColorscale colorscale;
	//protected StsColorscale currentColorscale = null;
	// these members are persistent, but not loaded from seis3d.txt.name file
	protected boolean readoutEnabled = false;

	transient public StsColorList seismicColorList = null;

	transient Attribute[] displayAttributes = null;
	transient Attribute currentAttribute = null;
	transient Attribute nullAttribute = new Attribute("none");

	// used by the ProgressTracker
//	transient protected double progressValue = 0.0;
//	transient protected String progressDescription = "";
//	transient protected StsProcessStepQueue processStepQueue = new StsProcessStepQueue();
	transient private StsActionListeners actionListeners = null;

	// the following are initialized after reading parameters file or reloading database
	transient protected StsSpectrumDialog spectrumDialog = null;
	transient protected boolean spectrumDisplayed = true;

//	static protected StsComboBoxFieldBean colorscalesBean;
//	static protected StsEditableColorscaleFieldBean colorscaleBean;

	// TODO move histogram inside colorscale or subclass a colorscale with histogram added.
	// TODO Since colorscaleBean is a static, instance data is added to beans when the panelObject is assigned,
	// TODO so histograms need to be handled this way as well.  Easiest thing is to fold them in with colorscale
	// TODO to which they belong anyways.  TJL 3/18/07
	abstract public void setDataHistogram();

    abstract public int getIntValue(int row, int col, int slice);
    abstract public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, boolean useByteCubes, float[] floatData);

    public StsSeismic()
	{
		super();
	}

	public StsSeismic(boolean persistent)
	{
		super(persistent);
    }

	/** volume has been removed from display; delete allocated planes, etc */
	public void deleteTransients()
	{
	}

	public boolean addToProject(boolean setAngle)
	{
		return currentModel.getProject().addToProject(this, setAngle);
	}

	public StsSpectrum getSpectrum()
	{
		return colorscale.getSpectrum();
	}

	public StsColorscale getColorscale()
	{
//		if(currentColorscale == colorscale) // colorscale is the seismic colorscale
		// would be nice to not have to set dataHistogram on every getColorscale().
		// presumably we are doing this because the colorscaleBean which has the histogram may not have been initialized
		// with histogram.  We should associate the histogram with the colorscale itself to avoid this confusion.
		// TJL 3/18/07
		setDataHistogram();
		return colorscale;
		// see comments at top of class about curentColorscale which is dysfunctional.
//        return currentColorscale;
	}

    public void initializeColorscaleActionListener()
    {
        StsColorscale colorscale = getColorscale();
		if (colorscale != null)
			colorscale.addActionListener(this);
    }

    public void reinitializeColorscale()
	{
		colorscale = null;
		initializeColorscale();
	}

	public boolean setGLColorList(GL gl, boolean nullsFilled, int dir, int shader)
	{
		if (dir == ZDIR && currentAttribute != null && currentAttribute != nullAttribute)
			return currentAttribute.setGLColorList(gl, nullsFilled, shader);
		else if (seismicColorList != null)
			return seismicColorList.setGLColorList(gl, nullsFilled, shader);
        else
			return false;
	}

	public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
		if(seismicColorList != null)
			return seismicColorList.setGLColorList(gl, nullsFilled, shader);
        else
			return false;
	}

	public int getColorDisplayListNum(GL gl, boolean nullsFilled)
	{
		if (!seismicColorList.setGLColorList(gl, nullsFilled, StsJOGLShader.NONE))
			return 0;
		return seismicColorList.colorListNum;
	}

    public void setColorscale(StsColorscale colorscale)
	{
		seismicColorList.setColorscale(colorscale);
		this.colorscale = colorscale;
		currentModel.win3dDisplayAll();
	}

	// new versions of get/set colorscale: not debugged yet.  TJL 9/24/06
/*
	  public void setColorscale(StsColorscale colorscale)
	  {
		StsColorList currentColorList = null;
		if(currentAttribute == nullAttribute)
	  currentColorList = seismicColorList;
else
	  currentColorList = currentAttribute.colorList;
		currentColorList.setColorscale(colorscale);
		this.colorscale = colorscale;
		currentModel.viewObjectChanged(this);
	  }

	  public StsColorscale getColorscale()
	  {
		if(currentAttribute == nullAttribute)
		{
	  setDataHistogram();
	  return seismicColorList.colorscale;

		}
else
	  return currentAttribute.colorList.colorscale;
	  }
*/
/*
	  public void itemStateChanged(ItemEvent e)
	  {
		if (e.getItem() instanceof StsColorscale)
		{
	  int id = e.getID();
	  if (id == StsColorscale.VOXELS_CHANGED)
	  {
		voxelVertexArraySets = null;
	  }
	 else
	  {
		currentModel.displayIfCursor3dObjectChanged(this);
		fireItemStateChanged(e);
	  }
		}
else
		{
	  ItemEvent event = null;
	  fireItemStateChanged(event);
		}
		return;
	  }
*/
	public void resetColors()
	{
		for (int i = 0; i < currentModel.viewPersistManager.families.length; i++)
		{
			StsWin3dBase[] windows = currentModel.getWindows(i);
			for (int n = 0; n < windows.length; n++)
			{
				StsWin3dBase window = (StsWin3dBase)windows[n];
				window.getCursor3d().objectChanged(this);
			}
		}
	}

	public void initializeColorscale()
	{
		try
		{
			if (colorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSeismicClass seismicClass = (StsSeismicClass)getCreateStsClass();
				colorscale = new StsColorscale("Seismic", spectrumClass.getSpectrum(seismicClass.getSeismicSpectrumName()), dataMin, dataMax);
				colorscale.setEditRange(dataMin, dataMax);
			}
            seismicColorList = new StsColorList(colorscale);
			colorscale.addActionListener(this);
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicVolume.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

    public void setDataRange(float dataMin, float dataMax)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        if(colorscale != null)
            colorscale.setRange(dataMin, dataMax);
        else
            initializeColorscale();
    }

    public StsColorscale getColorscaleWithName(String name)
	{
		if (name.equals(colorscale.getName()))
			return colorscale;
        else
			return null;
	}

	public void setSpectrumDialog(StsSpectrumDialog spectrumDialog)
	{
		this.spectrumDialog = spectrumDialog;
	}

	public boolean getReadoutEnabled()
	{
		return readoutEnabled;
	}

	public void setReadoutEnabled(boolean enabled)
	{
		readoutEnabled = enabled;
		return;
	}

	public final Color getColor(float[] xyz)
	{
		int index = getIntValue(xyz);
		return getColor(index);
	}

	public final Color getColor(int colorIndex)
	{
		return colorscale.getColor(colorIndex);
	}

	public final Color getColor(int row, int col, float z)
	{
		int index = getIntValue(row, col, z);
		return getColor(index);
	}

	public final Color getColor(float x, float y, float z)
	{
		int index = getIntValue(x, y, z);
		return getColor(index);
	}

	public StsColor getStsColor(float rowF, float colF, int plane)
	{
		int row = StsMath.roundOffInteger(rowF);
		int col = StsMath.roundOffInteger(colF);
		int v = getIntValue(row, col, plane);
		return colorscale.getStsColor(getIntValue(row, col, plane));
	}

	public float getValue(float x, float y, float z)
	{
		int row = getNearestBoundedRowCoor(y);
		int col = getNearestBoundedColCoor(x);
		int plane = getNearestBoundedSliceCoor(z);
		float f = (float)getIntValue(row, col, plane);
		return dataMin + (f / 254) * (dataMax - dataMin);
	}

	public float getValue(int row, int col, int slice)
	{
		float f = (float)getIntValue(row, col, slice);
		return dataMin + (f / 254) * (dataMax - dataMin);
	}

	public int getIntValue(float x, float y, float z)
	{
		int row = getNearestBoundedRowCoor(y);
		int col = getNearestBoundedColCoor(x);
		int plane = getNearestBoundedSliceCoor(z);
		return getIntValue(row, col, plane);
	}

	public int getIntValue(float[] xyz)
	{
		int row = getNearestBoundedRowCoor(xyz[1]);
		int col = getNearestBoundedColCoor(xyz[0]);
		int plane = getNearestBoundedSliceCoor(xyz[2]);
		return getIntValue(row, col, plane);
	}

	public int getIntValue(int row, int col, float z)
	{
		float kF = (z - zMin)/zInc;
		int k = (int)kF;
		float dk = kF - k;
		if (dk < 0.5f)
			return getIntValue(row, col, k);
		else
			return getIntValue(row, col, k + 1);
	}

	public void setIsVisible(boolean isVisible)
	{
		this.isVisible = isVisible;
        dbFieldChanged("isVisible", isVisible);
		currentModel.viewObjectRepaint(this, this);
	}

	public boolean getIsVisible()
	{
		return isVisible;
	}

	public synchronized void addActionListener(ActionListener listener)
	{
		if(actionListeners == null) actionListeners = new StsActionListeners();
		actionListeners.add(listener);
	}

	public synchronized void removeActionListener(ActionListener listener)
	{
		if(actionListeners == null) return;
		actionListeners.remove(listener);
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners == null) return;
		actionListeners.fireActionPerformed(e);
	}

	public class Attribute implements ActionListener
	{
		String name;
		double minValue;
		double maxValue;
		byte[] bytes;
		StsColorList colorList = null;
		StsTextureTiles textureTiles = null;
		boolean textureChanged = true;

		Attribute(String name)
		{
			this.name = name;
		}

		boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
		{
			if (colorList == null)
			{
				StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
				StsColorscale colorscale = new StsColorscale("Attribute", spectrum, 0.0f, 1.0f);
				colorList = new StsColorList(colorscale);
				colorscale.addActionListener(this);
			}
			colorList.setGLColorList(gl, nullsFilled, shader);
			return true;
		}

		public void actionPerformed(ActionEvent actionEvent)
		{
			colorList.setColorListChanged(true);
		}

		public String toString()
		{
			return name;
		}

		protected ByteBuffer getByteBuffer()
		{
			try
			{
				maxValue = -StsParameters.largeDouble;
				minValue = StsParameters.largeDouble;
				double[] values;
				values = getAttributeArray();
				int nValues = nRows * nCols;
				for (int n = 0; n < nValues; n++)
				{
					if (values[n] < minValue)
						minValue = values[n];
					else if (values[n] > maxValue)
						maxValue = values[n];
				}
				double scale = 254.0f / (maxValue - minValue);
				bytes = new byte[nValues];
				for (int n = 0; n < nValues; n++)
					bytes[n] = StsMath.unsignedIntToUnsignedByte254((int)(scale * (values[n] - minValue)));
				colorList.colorscale.setRange((float)minValue, (float)maxValue);
				return ByteBuffer.wrap(bytes);
			}
			catch (Exception e)
			{
				StsException.outputException("StsSeismicVolume.Attribute.getByteArray() failed for attribute: " + name,
													  e, StsException.WARNING);
				return null;
			}
		}

		private double[] getAttributeArray()
		{
			StsMappedDoubleBuffer attributeBuffer = null;

			try
			{
				attributeBuffer = StsMappedDoubleBuffer.constructor(stsDirectory, attributeFilename, "r");
				if (attributeBuffer == null)
					return null;
				int attributeIndex = getAttributeIndex(name);
				if (attributeIndex == -1)
				{
					new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find " + name + " attribute in attribute names list.");
					return null;
				}
				int nValues = nRows * nCols;
				long offset = (long)attributeIndex * nValues;
				if(!attributeBuffer.map(offset, nValues)) return null;
				double[] values = new double[nValues];
				attributeBuffer.get(values);
				return values;
			}
			catch (Exception e)
			{
				StsException.systemError("StsSeismicVolume.GetAttributeBuffer() failed.");
				return null;
			}
			finally
			{
				if (attributeBuffer != null)
					attributeBuffer.close();
			}
		}
	}
}
