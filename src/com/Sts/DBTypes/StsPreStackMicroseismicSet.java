package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Interpolation.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 8, 2006
 * Time: 8:30:45 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsPreStackMicroseismicSet extends StsPreStackLineSet implements StsTreeObjectI
{   
	static public final String SEISMIC_STRING = "Seismic";
	static public final String SEMBLANCE_STRING = "Semblance";
	
    static public final byte CVS = 0;
    static public final byte SEMBLANCE = 1;
    
    public String segyDirectory;
    static StsObjectPanel objectPanel = null;
    
	static public StsFieldBean[] displayFields = null;
	static public StsFieldBean[] propertyFields = null;
	
	public StsFilterProperties filterProperties = null;
	public StsAGCPreStackProperties agcProperties = null;
	public StsSemblanceDisplayProperties semblanceDisplayProperties = null;
	public StsSemblanceDisplayProperties semblanceEditProperties = null;
	public StsSemblanceComputeProperties semblanceComputeProperties = null;
	public StsSemblanceRangeProperties semblanceRangeProperties = null;
    public StsCVSProperties cvsProperties = null;
    
	transient public StsPreStackMicroseismicSetClass microseismicSetClass;
	
	protected StsObjectRefList colorscales = null;
	protected StsColorscale currentColorscale = null;
	transient public StsColorList seismicColorList = null;
	transient public StsColorList semblanceColorList = null;
	
	static protected StsEditableColorscaleFieldBean colorscaleBean;
	static protected StsComboBoxFieldBean colorscalesBean;
	
	transient protected StsColorscale seismicColorscale;
	transient protected StsColorscale semblanceColorscale;
	
	static GatherComparator gatherComparator = new GatherComparator();
	
	public StsMicroseismicGatherFile[] microGatherFiles = new StsMicroseismicGatherFile[0];
	public StsSEGYFormat segyFmt = null;
	
    public StsPreStackMicroseismicSet()
	{
    	;
    }

	public StsPreStackMicroseismicSet(boolean persistent)
	{
		super(persistent);
	}
	
	public void setMicroSegyFormat(StsSEGYFormat fmt)
	{
		segyFmt = fmt;
	}

	public boolean initialize(StsModel model)
	{
		try
		{
			microseismicSetClass = (StsPreStackMicroseismicSetClass) getStsClass();

            setupColorscales();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackMicroseismicSet.classInitialize(StsModel) failed.", e, StsException.WARNING);
			return false;
		}
	}
	
    public void initialize(float zMin, float zMax, float dataMin, float dataMax)
    {
    	setZMin(zMin);
    	setZMax(zMax);
    	setDataMin(dataMin);
    	setDataMax(dataMax);
        initialize();
	}
    
    public void initialize()
    {
        microseismicSetClass = (StsPreStackMicroseismicSetClass) getStsClass();
        
        filterProperties = new StsFilterProperties(this, microseismicSetClass.defaultFilterProperties, "filterProperties");
        semblanceDisplayProperties = new StsSemblanceDisplayProperties(this, microseismicSetClass.defaultSemblanceDisplayProperties, StsSemblanceDisplayProperties.DISPLAY_MODE, "semblanceDisplayProperties");
        semblanceEditProperties = new StsSemblanceDisplayProperties(this, microseismicSetClass.defaultSemblanceEditProperties, StsSemblanceDisplayProperties.EDIT_MODE, "semblanceEditProperties");
        agcProperties = new StsAGCPreStackProperties(this, microseismicSetClass.defaultAGCProperties, "agcProperties");
        cvsProperties = new StsCVSProperties(this, microseismicSetClass.defaultCVSProperties, "cvsProperties");
        semblanceRangeProperties = new StsSemblanceRangeProperties(this, currentModel, "semblanceRangeProperties");
        semblanceComputeProperties = new StsSemblanceComputeProperties(this, microseismicSetClass.defaultSemblanceComputeProperties, "semblanceComputeProperties");
		
        initializeColorscale();
		addToModel();
		semblanceRangeProperties.initializeSemblanceZRange(currentModel);        
		currentModel.setCurrentObject(this);
		initialize(currentModel);
	}
    
    public boolean addMicroseismicGather(String file, long start, long end)
    {
        try
        {
        	if(start < getBornDateLong() || getBornDateLong() == -1L)
        		setBornDate(start);
        	if(end > getDeathDateLong() || getDeathDateLong() == -1L)
        		setDeathDate(end);
        	
        	StsMicroseismicGatherFile gather = new StsMicroseismicGatherFile(file, start, end);
        	microGatherFiles = (StsMicroseismicGatherFile[]) StsMath.arrayAddSortedElement(microGatherFiles, gather, gatherComparator);        	
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPreStackMicroseismicSet.addMicroseismicGather() failed.", e, StsException.WARNING);
            return false;
        }
	} 
    
    public boolean isInTimeRange(long time)
    {
    	if(microGatherFiles == null) return false;
    	if((time > microGatherFiles[0].getStartTime()) && (time < microGatherFiles[microGatherFiles.length-1].getEndTime()))
    			return true;
    	return false;
    }
    
    public StsMicroseismicGatherFile getGather(long time)
    {
    	for(int i=0; i<microGatherFiles.length; i++)
        	if((time > microGatherFiles[i].getStartTime()) && (time < microGatherFiles[i].getEndTime()))
    			return microGatherFiles[i];
    	return null;
    }
    public void setupColorscales()
	{
		seismicColorscale = getColorscaleWithName(SEISMIC_STRING);
		if (seismicColorscale != null)
		{
			seismicColorList = new StsColorList(seismicColorscale);
			seismicColorscale.addActionListener(this);
		}
		semblanceColorscale = getColorscaleWithName(SEMBLANCE_STRING);
		if (semblanceColorscale != null)
		{
			semblanceColorList = new StsColorList(semblanceColorscale);
			semblanceColorscale.addActionListener(this);
		}
    } 
    
	public StsColorscale getColorscaleWithName(String name)
	{
		int nColorscales = colorscales.getSize();
		for(int n = 0; n < nColorscales; n++)
		{
			StsColorscale colorscale = (StsColorscale)colorscales.getElement(n);
			if(colorscale.getName().equals(name)) return colorscale;
		}
		return null;
	}
	
    public boolean addToProject()
	{
        return true;
	}

    public StsFieldBean[] getDisplayFields()
	{
		try
		{
			if (displayFields == null)
			{
				colorscalesBean = new StsComboBoxFieldBean(StsPreStackMicroseismicSet.class, "currentColorscale", "Colorscales");				
				colorscaleBean = new StsEditableColorscaleFieldBean(StsPreStackMicroseismicSet.class, "currentColorscale");
				
                StsIntFieldBean traceThresholdBean = new StsIntFieldBean(StsPreStackMicroseismicSet.class, "traceThreshold", 0, 100, "Minimum Traces Required:", true);
                traceThresholdBean.setToolTipText("Minimum number of traces required in a gather for velocity analysis");

                displayFields = new StsFieldBean[]
				{
					new StsBooleanFieldBean(StsPreStackMicroseismicSet.class, "isVisible", "Enable"),
					new StsBooleanFieldBean(StsPreStackMicroseismicSet.class, "readoutEnabled", "Mouse Readout"),
                    colorscalesBean,
					colorscaleBean,
					new StsButtonFieldBean("Semblance Display Properties", null, StsPreStackMicroseismicSet.class, "semblanceDisplayPropertiesDialog"),
					new StsButtonFieldBean("Gather Display Properties", null, StsPreStackMicroseismicSet.class, "gatherDisplayPropertiesDialog")
				};
			}
			colorscalesBean.setListItems(colorscales.getElements());
			colorscalesBean.setSelectedItem(currentColorscale);
			colorscaleBean.setValueObject(currentColorscale);
			colorscaleBean.setHistogram(dataHist);
			return displayFields;
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackMicroseismicSet.getDisplayFields() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public StsFieldBean[] getPropertyFields()
	{
		try
		{
			if (propertyFields == null)
			{
				propertyFields = new StsFieldBean[]
					{
					new StsStringFieldBean(StsPreStackMicroseismicSet.class, "name", true, "Name"),
					new StsFloatFieldBean(StsPreStackMicroseismicSet.class, "zMin", false, "Min Z or T"),
					new StsFloatFieldBean(StsPreStackMicroseismicSet.class, "zMax", false, "Max Z or T"),
					new StsFloatFieldBean(StsPreStackMicroseismicSet.class, "dataMin", false, "Data Min"),
					new StsFloatFieldBean(StsPreStackMicroseismicSet.class, "dataMax", false, "Data Max"),
					new StsFloatFieldBean(StsPreStackMicroseismicSet.class, "dataAvg", false, "Data Avg"),
		            new StsDateFieldBean(StsPreStackMicroseismicSet.class, "bornDate", false, "Born Date:"),
		            new StsDateFieldBean(StsPreStackMicroseismicSet.class, "deathDate", false, "Death Date:"),
					new StsButtonFieldBean("Semblance AGC & filter properties", null, StsPreStackMicroseismicSet.class, "semblanceAgcAndFilterPropertiesDialog"),
					new StsButtonFieldBean("Gather AGC & filter properties", null, StsPreStackMicroseismicSet.class, "gatherAgcAndFilterPropertiesDialog")
				};
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackMicroseismicSet.getPropertyFields() failed.", e, StsException.WARNING);
			return null;
		}
		return propertyFields;
	}

	public void initializeColorscale()
	{
		try
		{
			if (seismicColorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSpectrum seismicSpectrum = spectrumClass.getSpectrum(microseismicSetClass.getSeismicSpectrumName());
				seismicColorscale = new StsColorscale(SEISMIC_STRING, seismicSpectrum, -1.0f, 1.0f);
			}
			addColorscale(seismicColorscale);

			if (semblanceColorscale == null)
			{
				StsSpectrumClass spectrumClass = currentModel.getSpectrumClass();
				StsSpectrum semblanceSpectrum = spectrumClass.getSpectrum(microseismicSetClass.getSemblanceSpectrumName());
				semblanceColorscale = new StsColorscale(SEMBLANCE_STRING, semblanceSpectrum, 0.0f, 1.0f);
                semblanceColorscale.setEditMax(0.25f);
                semblanceColorscale.setCompressionMode(StsColorscale.COMPRESSED);
                semblanceColorscale.setTransparencyMode(false);
			}
			addColorscale(semblanceColorscale);
        }
		catch (Exception e)
		{
			StsException.outputException("StsPreStackMicroseismicSet.initializeColorscale() failed.", e, StsException.WARNING);
		}
	}

	public void addColorscale(StsColorscale colorscale)
	{
		if(colorscales == null) 
			colorscales = StsObjectRefList.constructor(2, 2, "colorscales", this);
		colorscales.add(colorscale);
		setCurrentColorscale(colorscale);
		if(colorscalesBean == null) return;
		colorscalesBean.setListItems(colorscales.getElements());
		colorscalesBean.setSelectedItem(colorscale);
	}

	public void setCurrentColorscale(StsColorscale colorscale)
	{
		if (currentColorscale == colorscale)
		{
			return;
		}
		currentColorscale = colorscale;
		if(isPersistent()) dbFieldChanged("currentColorscale", currentColorscale);
		if(colorscaleBean != null) colorscaleBean.setValueObject(colorscale);
		if(colorscalesBean != null) colorscalesBean.setSelectedItem(colorscale);
	}

	public StsColorscale getCurrentColorscale()
	{
		return currentColorscale;
	}

	public StsColorscale getSeismicColorscale()
	{
		return seismicColorscale;
	}

    public StsColorscale getSemblanceColorscale()
	{
		return semblanceColorscale;
	}

	public void setSeismicColorscale(StsColorscale colorscale)
	{
		this.seismicColorscale = colorscale;
		currentModel.viewObjectRepaint(this, colorscale);
	}

	public void setStacksColorscale(StsColorscale colorscale)
	{
		this.semblanceColorscale = colorscale;
		currentModel.viewObjectRepaint(this, colorscale);
	}
    public boolean anyDependencies()
    {
        return false;
    }
    
    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
			colorscaleBean.setHistogram(dataHist);
	}
	
    public void treeObjectSelected()
    {
        ;
    }
    
    public Object[] getChildren()
    {
        return new Object[0];
    }

	public void gatherPropertiesDialog()
	{
		gatherPropertiesDialog(currentModel.win3d);
	}

	public void gatherPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {getWiggleDisplayProperties(), agcProperties, filterProperties}, "Gather Properties", false);	}

	public void gatherDisplayPropertiesDialog()
	{
		gatherDisplayPropertiesDialog(currentModel.win3d);
	}

	public void gatherDisplayPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {getWiggleDisplayProperties()}, "Gather Display Properties", false);
	}

	public void gatherAgcAndFilterPropertiesDialog()
	{
		gatherAgcAndFilterPropertiesDialog(currentModel.win3d);
	}

	public void gatherAgcAndFilterPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {agcProperties, filterProperties}, "Gather AGC & Filter Properties", false);
	}
/*
	public void semblancePropertiesDialog()
	{
		semblancePropertiesDialog(currentModel.win3d);
	}

	public void semblancePropertiesDialog(Frame frame)
	{
		StsBatchOkApplyCancelDialog dialog = new StsBatchOkApplyCancelDialog(frame, this, this.SEMBLANCE, new StsDialogFace[]
								   {semblanceComputeProperties, semblanceRangeProperties, getSemblanceDisplayProperties(),
                                    agcProperties, filterProperties},
								   "Semblance Properties", false);
	}
*/
    //TODO implement this
    public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, boolean useByteCubes, float[] floatData)
    {
        return false;
    }

    public int getIntValue(int row, int col, int slice)
    {
        return 0;
    }
    
	public void semblanceDisplayPropertiesDialog()
	{
		semblanceDisplayPropertiesDialog(currentModel.win3d);
	}

	public void semblanceDisplayPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {semblanceComputeProperties, semblanceRangeProperties, getSemblanceDisplayProperties()}, "Semblance Display Properties", false);
	}

	public void semblanceAgcAndFilterPropertiesDialog()
	{
		semblanceAgcAndFilterPropertiesDialog(currentModel.win3d);
	}

	public void semblanceAgcAndFilterPropertiesDialog(Frame frame)
	{
		new StsOkApplyCancelDialog(frame, new StsDialogFace[]
								   {agcProperties, filterProperties}, "Semblance AGC and Filter Properties", false);
	}
/*
    public void cvsDisplayPropertiesDialog()
    {
        cvsPropertiesDialog(currentModel.win3d);
	}

	public void cvsPropertiesDialog(Frame frame)
	{
		StsBatchOkApplyCancelDialog dialog = new StsBatchOkApplyCancelDialog(frame, this, this.CVS, new StsDialogFace[]
		   {cvsProperties, semblanceRangeProperties, getSemblanceDisplayProperties(), getWiggleDisplayProperties(),
                   agcProperties, filterProperties}, "CVS Properties", false);
    }
*/
	public StsSemblanceDisplayProperties getSemblanceDisplayProperties()
	{
		return semblanceDisplayProperties;
	}

    public void displayWiggleProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {getWiggleDisplayProperties()}, "Wiggle Properties", false);
    }
    public void displaySemblanceProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {getSemblanceDisplayProperties() }, "Velocity Display Properties", false);
    }
    public void displayAGCProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {agcProperties}, "AGC Properties", false);
    }
    public void displayFilterProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {filterProperties}, "Filter Properties", false);
    }
    public void displayCVSProperties()
    {
        new StsOkApplyCancelDialog(currentModel.win3d, new StsDialogFace[] {cvsProperties}, "Edit CVS/VVS Properties", false);
    }
	public StsSemblanceComputeProperties getSemblanceComputeProperties()
	{
		return semblanceComputeProperties;
	}
    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return microseismicSetClass.getWiggleDisplayProperties();
    }
    
    public void updatePreStackVolumes(StsRadialInterpolation interpolation)  {  }
    public int getVolumeRowColIndex(int row, int col) {	return 1; }
    public String getEmptyGatherDescription(int row, int col) { return new String("Microseismic Gather"); }
    public String getFullGatherDescription(int row, int col, int nTraces) { return new String("Microseismic Gather with " + nTraces + " traces.");}
	public Component getPanel(boolean val) { return getPanel(); }
	public Component getPanel()
	{
		 exportPanel = new StsSeismicExportPanel2d(currentModel, this, "Export Microseismic Gather", true);
		 return exportPanel;
	}    
    public int[] adjustLimitRowCol(int nRow, int nCol)
    {
        nRow = StsMath.minMax(nRow, 0, nRows-1);
        nCol = StsMath.minMax(nCol, 0, lines[nRow].nCols-1);
        return new int[] { nRow, nCol };
    }
    public Object[] getAvailableModelsAndVolumesList(){ return null; }
    public Object[] getAvailableModelsList()
	{
		Object[] models = currentModel.getTrimmedList(StsPreStackVelocityModel.class);
		if (models.length == 0)
		{
			models = new Object[] {NO_MODEL};
		}
		return models;
	}
	public void setInputVelocityModel(Object object)
	{
		if (object instanceof StsPreStackVelocityModel2d)
		{
			dbFieldChanged("inputVelocityModel", inputVelocityModel, object);
		}
		else
		{
			dbFieldChanged("inputVelocityModel", inputVelocityModel, null);
		}
	}
	public void initializeVelStatusPoints()
	{
        ;
	}	
    public boolean allocateVolumes(String mode, boolean loadFromFile)
	{
    	return true;
    }
	public boolean isPlaneOK(byte displayMode, int dir, int nPlane)
	{
		return false;
	}
	public void computePreStackVolume(byte displayMode)
	{
	}	
    public StsPreStackVelocityModel constructVelocityModel()
    {
        return null;
    }

    public int[] getRowColFromCoors(float x, float y)
    {
        float dSq = StsParameters.largeFloat;
        int nearestRow = -1;
        int nearestCol = -1;
        for(int row = 0; row < lines.length; row++)
        {
            StsPreStackLine line = lines[row];
            int nCols = line.nCols;
            float[] cdpX = line.cdpX;
            float[] cdpY = line.cdpY;
            for (int col = 0; col < nCols; col++)
			{
				float dx = (cdpX[col] - x);
                float dy = (cdpY[col] - y);
                float dSqNew = dx*dx + dy*dy;
                if(dSqNew < dSq)
                {
                    nearestRow = row;
                    nearestCol = col;
                    dSq = dSqNew;
                }
            }
        }
        return new int[] { nearestRow, nearestCol } ;
    }

	/** for prestack 2d, each line is a row and col is gather position along the line  from 0 to nGathers-1 */
	public StsPreStackLine getDataLine(int row, int col)
	{
		if (lines == null)
		{
			return null;
		}

		int nLine = row;
		int nGather = col;

		int nLines = lines.length;
		if (nLine < 0 || nLine >= nLines)
		{
			return null;
		}

		StsPreStackLine line = lines[nLine];
		int compare = line.compareTo(nGather);
		if (compare == 0)
			return line;
		else
			return null;
	}
    public int getNColsForRow(int row)
    {
        if(row >= nRows) return 0;
        return lines[row].nCols;
    }
/*
    public StsMicroseismicSuperGather getGather()
	{
		return getSuperGather(currentModel.glPanel3d);
	}
    
	public StsMicroseismicSuperGather getSuperGather(StsGLPanel3d glPanel3d)
	{
		return getSuperGather(glPanel3d, null);
	}    
	
	public StsMicroseismicSuperGather getSuperGather(StsGLPanel3d glPanel3d, StsMicroseismicGatherFile gather)
	{
        if(glPanel3d == null) return null;
        if(gather == null) 
        	gather = microGatherFiles[0];
        
        StsMicroseismicSuperGather microGather = gather.buildSuperGather(this);
		return microGather;
	}
	*/
    static class GatherComparator implements Comparator
	{
		GatherComparator()
		{
		}

		public int compare(Object o1, Object o2)
		{
			//TJL long z1 = ((StsMicroseismicGather)o1).getStartTime();
			//TJL long z2 = ((StsMicroseismicGather)o2).getStartTime();
			//TJL if(z2 > z1) return 1;
			//TJL else if(z2 < z1) return -1;
			//TJL else
            return 0;
		}
		public boolean equals(Object o1)
		{
			return false;
		}
	}   
}