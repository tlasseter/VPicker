package com.Sts.DBTypes;

import com.Sts.Interfaces.StsFractureDisplayable;
import com.Sts.Interfaces.StsTimeEnabled;
import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.Interfaces.StsView3dDisplayable;
import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.Types.StsBoundingBox;
import com.Sts.Types.StsGolderFracture;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.StsObjectPanel;
import com.Sts.UI.StsSelectable;
import com.Sts.Utilities.*;

import javax.media.opengl.GL;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 5:45:18 PM
 * To change this tempfclate use File | Settings | File Templates.
 */
public class StsGolderFractureSet extends StsBoundingBox implements StsSelectable, StsTreeObjectI, StsTimeEnabled, StsView3dDisplayable
{
	public int setNumber;
    protected StsGolderFracture[] fractures;
    protected int nFractures = 0;
    protected String[] propertyNames;
	protected float[][] propertyRanges;
    protected byte zDomainSupported = StsParameters.TD_TIME_DEPTH;
    protected float colorMin = StsParameters.largeFloat;
    protected float colorMax = -StsParameters.largeFloat;
    protected String spectrumName = "Basic";
    protected StsColorscale colorscale;
    protected boolean displaySectionEdges = false;
   /** Index of color or vector to color events: 0=NONE, 1=STAGE, 2-33 = 32 colors, 34+ = valueVectorIndex - 33 */
    protected int colorItemIndex = COLOR_NONE;
	int fracsDisplayed = 0;

	/** colorByList consists of colors based on states (NONE, STAGE, etc), properties (AMPLITUDE, S/N, etc), and general colors */
    transient String[] colorByList;
	transient StsColor[] fractureColors;

    transient protected byte zDomainDisplayed = StsParameters.TD_NONE;
    transient boolean useDisplayLists;
    transient boolean usingDisplayLists = true;

    transient boolean colorscaleChanged = false;
	/** indicates colorVector has been changed */
    transient protected boolean colorChanged = true;

    static protected StsEditableColorscaleFieldBean colorscaleBean;
    static protected StsComboBoxFieldBean colorByBean;
    static protected StsFloatFieldBean colorMinBean;
    static protected StsFloatFieldBean colorMaxBean;
    static StsDateFieldBean bornField = null;
    static StsDateFieldBean deathField = null;

	static public final int COLOR_NONE = 0, COLOR_SEQUENCE = 1;
    static public final String COLOR_NONE_STRING = "None", COLOR_SEQUENCE_STRING = "Sequence";
    static public String[] colorActionsList;
    static { colorActionsList = new String[] { COLOR_NONE_STRING, COLOR_SEQUENCE_STRING };  }

	static StsGolderFractureSetClass fractureSetClass = null;

    static final int fractureSetSizeIncrement = 1000;


    public StsGolderFractureSet()
    {
    }

    public StsGolderFractureSet(boolean persistent)
    {
        super(persistent);
    }

    public StsGolderFractureSet(String fractureSetName, String[] propertyNames, double xOrigin, double yOrigin)
    {
        super(false);
		setName(fractureSetName);
		this.propertyNames = propertyNames;
		setOrigin(xOrigin, yOrigin);
    }

    static public StsGolderFractureSet constructor()
    {
        try
        {
            StsGolderFractureSet fractureSet = new StsGolderFractureSet(false);
            fractureSet.initialize();
            return fractureSet;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsFractureSet.class, "constructor", e);
            return null;
        }
    }

	static public StsGolderFractureSet constructor(String fractureSetName, String[] propertyNames, double xOrigin, double yOrigin)
	{
       try
        {
			StsGolderFractureSet fractureSet = new StsGolderFractureSet(fractureSetName, propertyNames, xOrigin, yOrigin);
			fractureSet.initialize();
            return fractureSet;
		}
        catch (Exception e)
        {
            StsException.outputWarningException(StsFractureSet.class, "constructor(fractureSetName...)", e);
            return null;
        }
	}

    public boolean initialize(StsModel model)
    {
        checkCreateFractureSetClass();
        return true;
    }

    public void initialize()
    {
        checkCreateFractureSetClass();
        StsColor color = new StsColor(fractureSetClass.getNextColor());
        setStsColor(color);
    }

	private void checkInitializeColorscale()
	{
		if (colorscale == null)
		{
			spectrumName = getGolderFractureSetClass().getDefaultSpectrumName();
			colorscale = new StsColorscale("GolderFractureSet", currentModel.getSpectrum(spectrumName), 0.0f, 1.0f);
		}
		colorscale.addActionListener(this);
	}

    public int getNFractures()
    {
        return nFractures;
    }

    public String[] getPropertyNames() { return propertyNames; }
    public void setPropertyNames(String[] names){ propertyNames = names; }

	    public String[] getColorByList()
    {
        if(colorByList != null) return colorByList;
        return initializeColorByList();
    }

    public String[] initializeColorByList()
    {
        colorByList = (String[])StsMath.arrayAddArray(colorByList, colorActionsList);
        colorByList = (String[])StsMath.arrayAddArray(colorByList, propertyNames);
        colorByList = (String[])StsMath.arrayAddArray(colorByList, StsColor.colorNames32);
        return colorByList;
    }

    public void setColorBy(String colorBy)
    {
        String[] colorByList = getColorByList();
		if(StsStringUtils.stringsEqual(colorBy, colorByList[colorItemIndex])) return;
        for(int i = 0; i < colorByList.length; i++)
        {
            if(StsStringUtils.stringsEqual(colorBy, colorByList[i]))
            {
                if(colorItemIndex == i) return;
                colorItemIndex = i;
                colorChanged = true;
                dbFieldChanged("colorItemIndex", colorItemIndex);
                currentModel.viewObjectRepaint(this, this);
                return;
            }
        }
        colorItemIndex = 0;
    }

    public String getColorBy()
    {
        String[] colorByList = getColorByList();
        return colorByList[colorItemIndex];
    }
    private void setColors()
    {
		int nColorActions = colorActionsList.length;
		int nValues = fractures.length;
		int nProperties = propertyNames.length;
        fractureColors = new StsColor[fractures.length];
 		if(colorItemIndex == COLOR_NONE)
		{
			for(int n = 0; n < nValues; n++)
				fractureColors[n] = StsColor.GREY;
		}
        else if(colorItemIndex == COLOR_SEQUENCE)
		{
			StsColor color;
			int objectIndex = getIndex();
			if(objectIndex != -1)
				color = StsColor.colors32[objectIndex % 32];
			else
				color = StsColor.GREY;
			for(int n = 0; n < nValues; n++)
				fractureColors[n] = color;
		}
		else if(colorItemIndex >= colorActionsList.length + nProperties)
		{
			StsColor color = StsColor.colors32[(colorItemIndex - colorActionsList.length - nProperties)%32];
			for(int n = 0; n < nValues; n++)
				fractureColors[n] = color;
		}
        else // use a property
        {
			int nProperty = colorItemIndex - nColorActions;
			float[] propertyRange = propertyRanges[nProperty];
			float scale = StsMath.floatToUnsignedByteScale(propertyRange);
			float offset = StsMath.floatToUnsignedByteScaleOffset(scale, propertyRange[0]);
            for(int n = 0; n < nFractures; n++)
			{
				float value = fractures[n].getPropertyValue(nProperty);
				int index = StsMath.floatToInt254WithScale(value, scale, offset);
				fractureColors[n] = colorscale.getStsColor(index);
        	}
		}
    }

    public int getFracsDisplayed()
    {
        return fracsDisplayed;
    }

    public void setFracsDisplayed(int num)
    {
        fracsDisplayed = num;
    }
    public int getFracCount()
    {
        return nFractures;
    }
    public void setFracCount(int num) {};

    boolean decimationChanged = false;

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;
    static protected StsObjectPanel objectPanel = null;
    static protected StsIntFieldBean totalFracsBean;
    static protected StsIntFieldBean fracsDisplayedBean;
    static protected StsIntFieldBean decimationBean;
    static protected StsFloatFieldBean verticalScaleBean;
    transient private int listNum = 0; // display list number (>0)
    private int decimationPercent = 0;
    private float verticalScale = 1.0f;

    public int getDecimationPercent()
    {
        return decimationPercent;
    }

    public void setDecimationPercent(int percent)
    {
        if (this.decimationPercent != percent)
        {
            decimationChanged = true;
            decimationPercent = percent;
            int numToDisplay = ((100 - decimationPercent) * nFractures) / 100;
            if (fracsDisplayedBean != null)
                fracsDisplayedBean.setValue(numToDisplay);
            dbFieldChanged("decimationPercent", decimationPercent);
            currentModel.viewObjectRepaint(this, this);
        }
    }
    public float getVerticalScale()
    {
        return verticalScale;
    }

    public void setVerticalScale(float scale)
    {
        if (this.verticalScale != scale)
        {
            verticalScale = scale;
            decimationChanged = true;
            dbFieldChanged("verticalScale", verticalScale);
            currentModel.viewObjectRepaint(this, this);
        }
    }

    /**
     * Get the fracture set class containing this object
     *
     * @return
     */
    static public StsGolderFractureSetClass getGolderFractureSetClass()
    {
        return (StsGolderFractureSetClass) currentModel.getCreateStsClass(StsGolderFractureSet.class);
    }

    public void display(StsGLPanel glPanel)
    {
        if (nFractures == 0) return;
        if(getGolderFractureSetClass().getEnableTime() && !isAlive(currentModel.getProject().getProjectTime()))
            return;

		GL gl = glPanel.getGL();

        if (currentModel.getProject().isTimeRunning())
            deleteDisplayLists(gl);

        if (!isVisible)
        {
            deleteDisplayLists(gl);
            return;
        }

        checkInitializeColorscale();

        if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
        {
            deleteDisplayLists(gl);
            return;
        }
        byte projectZDomain = currentModel.getProject().getZDomain();
        if (projectZDomain != zDomainDisplayed)
        {
            deleteDisplayLists(gl);
            zDomainDisplayed = projectZDomain;
        }
		if(colorscaleChanged || colorChanged)
		{
			setColors();
		}
        useDisplayLists = currentModel.useDisplayLists;
        if (!useDisplayLists && usingDisplayLists)
        {
            deleteDisplayLists(gl);
            usingDisplayLists = false;
        }
        if ((decimationChanged) || colorscaleChanged || colorChanged)
        {
            deleteDisplayLists(gl);
            decimationChanged = false;
            colorscaleChanged = false;
			colorChanged = false;
        }

        displayFractures((StsGLPanel3d) glPanel);
//        if (fracsDisplayedBean != null)
//        	fracsDisplayedBean.setValue(fracsDisplayed);
// 
    }

    private void deleteDisplayLists(GL gl)
    {
        if (listNum != 0)
        {
            gl.glDeleteLists(listNum, 1);
            listNum = 0;
        }
    }

    private void displayFractures(StsGLPanel3d glPanel)
    {
        //	boolean displaySectionEdges = false;
        GL gl = glPanel.getGL();

		double viewShift = getIndex();
        glPanel.setViewShift(gl, viewShift);

        if (useDisplayLists)
        {
            drawDisplayList(glPanel, displaySectionEdges);
        }
        else
        {
            if (listNum != 0) // delete existing display list
            {
                gl.glDeleteLists(listNum, 1);
                listNum = 0;
            }
            drawFractures(glPanel, displaySectionEdges);
        }
        glPanel.resetViewShift(gl);
    }

    private void drawDisplayList(StsGLPanel3d glPanel, boolean displaySectionEdges)
    {
        GL gl = glPanel.getGL();
        if (listNum == 0) // build display list
        {
            listNum = gl.glGenLists(1);
            if (listNum == 0)
            {
                StsMessageFiles.logMessage(
                    "System Error in StsGrid.displaySurface: " +
                        "Failed to allocate a display list");
                return;
            }

            gl.glNewList(listNum, GL.GL_COMPILE_AND_EXECUTE);
            drawFractures(glPanel, displaySectionEdges);
            gl.glEndList();
        }

        gl.glCallList(listNum);
    }

    private void drawFractures(StsGLPanel3d glPanel, boolean displaySectionEdges)
    {
        int numToDisplay = ((100 - decimationPercent) * nFractures) / 100;
        double step = (double) nFractures / numToDisplay;
        int count = 0;
        int displayed = 0;
        GL gl = glPanel.getGL();

        for(int i=0; i<fractures.length; i++)
        {
            int index = (int) (step / 2 + Math.floor(step * count));
            if(fractures[i].isHighlighted() || ((index == i) && (numToDisplay > 0)))
            {
                StsColor stsColor = getFractureColor(i);
                if(checkTime(fractures[i].getTime()))
                    fractures[i].display(glPanel, displaySectionEdges, stsColor);
                if(index == i)
                {
                    ++count;
                    displayed++;
                    if(count == numToDisplay) break;
                }
            }
        }
        setFracsDisplayed(displayed);
    }

	private StsColor getFractureColor(int index)
	{
        // Are the colors being overridden with active and inactive colors
        if(getGolderFractureSetClass().getColorActive())
        {
            if(fractures[index].highlight)
                return getGolderFractureSetClass().getActiveColor();
            else
                return getGolderFractureSetClass().getInactiveColor();
        }
        else
        {
            if(fractures[index].highlight)
		        return fractures[index].overrideFractureColor(getGolderFractureSetClass().getActiveColor());
            else
            {
                if(fractureColors == null)
			        return getGolderFractureSetClass().getInactiveColor();
		        else
			        return fractureColors[index];
            }
        }
	}

    /**
     * Determine if an event contained in this sensor meets the time criteria to display
     *
     * @param time - time of current fracture
     * @return success
     */
    private boolean checkTime(long time)
    {
        if (time == 0L) return true;

        if ((time > currentModel.getProject().getProjectTime())
            || (time == StsParameters.nullLongValue))
            return false;

        long duration = currentModel.getProject().getProjectTimeDuration();
        long projectTime = currentModel.getProject().getProjectTime();

        if (duration != 0)
        {
            long endTime = projectTime - duration;
            if (time < endTime)
                return false;
        }
        return true;
    }

    public void checkCreateFractureSetClass()
    {
        if (fractureSetClass == null)
            fractureSetClass = (StsGolderFractureSetClass) currentModel.getCreateStsClass(StsGolderFractureSet.class);
    }

    public boolean delete()
    {
        if (!super.delete()) return false;
        if (fractures != null) fractures = null;
        return true;
    }

    public void addFracture(StsGolderFracture fracture)
    {
        fractures = (StsGolderFracture[])StsMath.arrayAddElement(fractures, fracture, nFractures, fractureSetSizeIncrement);
        nFractures++;
		float[][] vertices = fracture.vertices;
        for (float[] vertex : vertices)
            addPoint(vertex, 0, 0);
    }

  public void addToModel()
    {
        currentModel.checkAddToCursor3d(this);
        super.addToModel();
    }


    public void setIsVisible(boolean value)
    {
        isVisible = value;
        dbFieldChanged("isVisible", isVisible);
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getIsVisible()
    {
        return isVisible;
    }

    /** Handle colorscale changes */
    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            colorscaleChanged = true;
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
    }

    /** Set the colorscale to be applied to the colorby curve */
    public void setColorscale(StsColorscale colorscale)
	{
		this.colorscale = colorscale;
		colorscale.addActionListener(this);
	}

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            bornField = new StsDateFieldBean(StsGolderFractureSet.class, "bornDate", "Born Date:");
            deathField = new StsDateFieldBean(StsGolderFractureSet.class, "deathDate", "Death Date:");
            displayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsGolderFractureSet.class, "isVisible", "Enable"),
                    bornField,
                    deathField,
                    //new StsBooleanFieldBean(StsGolderFractureSet.class, "drawEdges", "Edges"),
                    totalFracsBean = new StsIntFieldBean(StsGolderFractureSet.class, "fracCount", "Number Fractures:"),
                    fracsDisplayedBean = new StsIntFieldBean(StsGolderFractureSet.class, "fracsDisplayed", "Enable Fractures:"),
                    new StsBooleanFieldBean(StsGolderFractureSet.class, "displaySectionEdges", "Display Edges"),
                    decimationBean = new StsIntFieldBean(StsGolderFractureSet.class, "decimationPercent", 0, 100, "% Decimation::", true),
                    verticalScaleBean = new StsFloatFieldBean(StsGolderFractureSet.class, "verticalScale", true, "Vertical Scale",true),
                    colorByBean = new StsComboBoxFieldBean(StsGolderFractureSet.class, "colorBy", "Color By:", "colorByList"),
                    colorMinBean = new StsFloatFieldBean(StsGolderFractureSet.class, "colorMinimum", true, "Color Min"),
                    colorMaxBean = new StsFloatFieldBean(StsGolderFractureSet.class, "colorMaximum", true, "Color Max"),
                    colorscaleBean = new StsEditableColorscaleFieldBean(StsGolderFractureSet.class, "colorscale")
                };
            fracsDisplayedBean.setEditable(false);
            totalFracsBean.setEditable(false);
            verticalScaleBean.setValueAndRangeFixStep(1.0f, 0.1f, 10.0f, 0.1f);
            //decimationBean.setEditable(false);
        }
        //fracsDisplayedBean.setValue(fracsDisplayed);

        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
//        if(propertyFields == null)
//        {
//            propertyFields = new StsFieldBean[]
//            {
//                new StsFloatFieldBean(this, "stimulatedRadius", 0.0f, 10000.0f, "Stimulated radius")
//            };
//        }
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }


    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel != null)
            return objectPanel;
        objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsGolderFractureSet").selected(this);
    }

    public String getLabel()
    {
        return new String("fracture set: " + getName());
    }

    public boolean getDisplaySectionEdges(){ return displaySectionEdges; }

    public void setDisplaySectionEdges(boolean displayEdges)
    {
        displaySectionEdges = displayEdges;
        dbFieldChanged("displaySectionEdges", displaySectionEdges);
        colorscaleChanged = true;        
        currentModel.viewObjectRepaint(this, this);
    }

    public float getColorMinimum(){ return colorMin; }

    public void setColorMinimum(float min)
    {
        colorMin = min;
        colorscale.setRange((float) colorMin, (float) colorMax);
        colorscaleChanged = true;
        dbFieldChanged("colorMin", colorMin);
        currentModel.viewObjectRepaint(this, this);
    }

    public void setColorMaximum(float max)
    {
        colorMax = max;
        colorscale.setRange((float) colorMin, (float) colorMax);
        colorscaleChanged = true;
        dbFieldChanged("colorMax", colorMax);
        currentModel.viewObjectRepaint(this, this);
    }

    public float getColorMaximum(){ return colorMax; }

	public StsColorscale getColorscale() { return colorscale; }

	public void finish()
	{
		fractures = (StsGolderFracture[])StsMath.trimArray(fractures, nFractures);
		computePropertyRanges();
	}
    
	/** There are nProperties number of properties.  Scan all fractures and set property ranges. */
	private void computePropertyRanges()
	{
		int nProperties = propertyNames.length;
		propertyRanges = new float[nProperties][2];
		for(int n = 0; n < nProperties; n++)
		{
			float firstPropertyValue = fractures[0].properties[n];
			propertyRanges[n][0] = firstPropertyValue;
			propertyRanges[n][1] = firstPropertyValue;
		}
		for(int i = 1; i < fractures.length; i++)
		{
			for(int n = 0; n < nProperties; n++)
			{
				float[] propertyRange = propertyRanges[n];
				float property = fractures[i].properties[n];
				if(property < propertyRange[0]) propertyRange[0] = property;
				else if(property > propertyRange[1]) propertyRange[1] = property;
			}
		}
	}

    public void setBornDate(String born)
    {
        super.setBornDate(born);
        bornField.setValueFromPanelObject(this);
        currentModel.viewObjectRepaint(this, this);
        return;
    }

    public void setDeathDate(String death)
    {
        super.setDeathDate(death);
        deathField.setValueFromPanelObject(this);
        currentModel.viewObjectRepaint(this, this);
        return;
    }

	public void highlightFracturesInsideLimits(StsDynamicSensor sensor, int index, float maxDistance, boolean azimuthLimit, float minAzimuth, float maxAzimuth)
	{
		if(!isInsideXYZ(sensor.getXYZ(index))) return;
		for(StsGolderFracture fracture : fractures)
			fracture.highlightIfInsideLimits(sensor, index, maxDistance, azimuthLimit, minAzimuth, maxAzimuth);
	}


	public boolean highlightIntersectedFractures(ArrayList<StsFractureDisplayable> interpretedFractures)
	{
		boolean intersected = false;
		int nInterpretedFractures = interpretedFractures.size();
		StsBoundingBox[] fractureBoundingBoxes = new StsBoundingBox[nInterpretedFractures];
		for(int n = 0; n < nInterpretedFractures; n++)
			fractureBoundingBoxes[n] = interpretedFractures.get(n).getBoundingBox();
		for(StsGolderFracture fracture : fractures)
		{
			for(StsFractureDisplayable interpretedFracture : interpretedFractures)
			{
				if(!interpretedFracture.getIsVisible()) continue;
				if(!interpretedFracture.intersects(fracture)) continue;
                fracture.setHighlight(interpretedFracture, 0);
				intersected = true;
			}
		}
		return intersected;
	}
}
