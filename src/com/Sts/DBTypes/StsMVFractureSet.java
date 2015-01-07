package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 5:45:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsMVFractureSet extends StsBoundingBox implements StsSelectable, StsTreeObjectI, StsTimeEnabled
{
    /**
     *
     */
    private static final long serialVersionUID = -474248455773803025L;
    protected StsMVFracture[] fractures = null;
    protected int nFractures = 0;
    protected String[] attributeNames = null;
    transient String[] attributeList = new String[]{StsMVFracture.S_NONE, StsMVFracture.S_APERTURE,
        StsMVFracture.S_LENGTH, StsMVFracture.S_DIP, StsMVFracture.S_AZIMUTH, StsMVFracture.S_ASPECTRATIO};
    protected String colorAttribute = StsMVFracture.S_NONE;
    protected int colorAttributeIdx = StsMVFracture.NONE;

    protected byte zDomainSupported = StsParameters.TD_TIME_DEPTH;
    transient protected byte zDomainDisplayed = StsParameters.TD_NONE;
    transient boolean useDisplayLists;
    transient boolean usingDisplayLists = true;
    StsColor stsColor;
    /** color of this section (and sectionPatch) */

    protected float colorMin = StsParameters.largeFloat;
    protected float colorMax = -StsParameters.largeFloat;
    protected String spectrumName = "Basic";
    protected StsColorscale colorscale;
    protected boolean displaySectionEdges = false;
    transient boolean colorscaleChanged = false;

    static protected StsEditableColorscaleFieldBean colorscaleBean;
    static protected StsComboBoxFieldBean colorByBean;
    static protected StsFloatFieldBean colorMinBean;
    static protected StsFloatFieldBean colorMaxBean;

    static final int fractureSetSizeIncrement = 1000;

    public int getNFractures()
    {
        return nFractures;
    }

    /** Get list of attributes */
    public String[] getAttributeNames()
    { return attributeNames; }

    public void setAttributeNames(String[] names){ attributeNames = names; }

    public String[] getAttributeList(){ return attributeList; }

    public String getColorBy(){ return colorAttribute; }

    public void setColorBy(String att)
    {
        colorAttribute = att;
        if (!colorAttribute.equalsIgnoreCase(StsMVFracture.S_NONE))
        {
            for (int i = 0; i < attributeList.length; i++)
            {
                if (attributeList[i].equalsIgnoreCase(colorAttribute))
                    colorAttributeIdx = i;
            }
            colorMin = StsParameters.largeFloat;
            colorMax = -StsParameters.largeFloat;
            for (int i = 0; i < nFractures; i++)
            {
                float val = fractures[i].getAttributeValue(colorAttributeIdx);
                if (val < colorMin)
                    colorMin = val;
                if (val > colorMax)
                    colorMax = val;
            }
            if (colorscale != null)
            {
                colorscale.setRange(colorMin, colorMax);
            }
            colorMinBean.setValueObject(colorMin);
            colorMaxBean.setValueObject(colorMax);
        }
        colorscaleChanged = true;
        currentModel.viewObjectRepaint(this, this);
    }

    int fracsDisplayed = 0;

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
            currentModel.viewObjectRepaint(this, this);
        }
    }

    private void initializeColorscale()
    {
        if (colorscale == null)
        {
            spectrumName = getMVFractureSetClass().getDefaultSpectrumName();
            colorscale = new StsColorscale("MVFractureSet", currentModel.getSpectrum(spectrumName), (float) colorMin, (float) colorMax);
        }
        colorscale.addActionListener(this);
    }

    /**
     * Get the fracture set class containing this object
     *
     * @return
     */
    static public StsMVFractureSetClass getMVFractureSetClass()
    {
        return (StsMVFractureSetClass) currentModel.getCreateStsClass(StsMVFractureSet.class);
    }

    /**
     * Compute the color for the supplied value
     *
     * @param value - attribute value
     * @return StsColor
     */
    private StsColor defineColor(float value)
    {
        // Define the color
        int bin = 0;
        double interval = (colorMax - colorMin) / (double) (colorscale.getNColors() - 1);
        if (!colorAttribute.equalsIgnoreCase(StsMVFracture.S_NONE))
        {
            bin = (int) (Math.abs((value - colorMin)) / interval);
            int colorIdx = bin % colorscale.getNColors();
            return colorscale.getStsColor(colorIdx);
        }
        else
            return stsColor;
    }

    static StsMVFractureSetClass fractureSetClass = null;

    public StsMVFractureSet()
    {
    }

    public StsMVFractureSet(boolean persistent)
    {
        super(persistent);
    }

    public void initialize()
    {
        checkCreateFractureSetClass();
        StsColor color = new StsColor(fractureSetClass.getNextColor());
        setStsColor(color);
        initializeAttributeList();
    }

    static public StsMVFractureSet constructor()
    {
        try
        {
            StsMVFractureSet fractureSet = new StsMVFractureSet(false);
            fractureSet.initialize();
            return fractureSet;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsFractureSet.class, "constructor", e);
            return null;
        }
    }

    public boolean initialize(StsModel model)
    {
        checkCreateFractureSetClass();
        initializeAttributeList();
        return true;
    }

    public void initializeAttributeList()
    {
        if (fractures != null)
            attributeNames = fractures[0].getAttributeNames();

        attributeList = (String[]) StsMath.arrayAddArray(attributeList, attributeNames);
    }

    public void finish()
    {
        fractures = (StsMVFracture[]) StsMath.trimArray(fractures, nFractures);
        addToModel();
    }

    public void display(StsGLPanel glPanel)
    {
        if (nFractures == 0) return;

        if (currentModel.getProject().isTimeRunning())
            listNum = 0;

        GL gl = glPanel.getGL();

        if (!isVisible)
        {
            deleteDisplayLists(gl);
            return;
        }

        if (colorscale == null) initializeColorscale();

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

        useDisplayLists = currentModel.useDisplayLists;
        if (!useDisplayLists && usingDisplayLists)
        {
            deleteDisplayLists(gl);
            usingDisplayLists = false;
        }
        if ((decimationChanged) || (colorscaleChanged))
        {
            deleteDisplayLists(gl);
            decimationChanged = false;
            colorscaleChanged = false;
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

        if (stsColor != null)
            stsColor.setGLColor(gl);
        else
            StsColor.PURPLE.setGLColor(gl);

        glPanel.setViewShift(gl, StsGraphicParameters.sectionShift);

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

        while (count < numToDisplay)
        {
            int index = (int) (step / 2 + Math.floor(step * count));
            if (!colorAttribute.equalsIgnoreCase(StsMVFracture.S_NONE))
            {
                stsColor = defineColor(fractures[index].getAttributeValue(colorAttributeIdx));
                //color.setGLColor(gl);
            }
            if (checkTime(fractures[index].getTime()))
                fractures[index].display(glPanel, displaySectionEdges, stsColor);
            ++count;
            displayed++;
        }
        setFracsDisplayed(displayed);
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
            fractureSetClass = (StsMVFractureSetClass) currentModel.getCreateStsClass(StsMVFractureSet.class);
    }


    public boolean delete()
    {
        if (!super.delete()) return false;
        if (fractures != null) fractures = null;
        return true;
    }

    public void addFracture(StsMVFracture fracture)
    {
        fractures = (StsMVFracture[])StsMath.arrayAddElement(fractures, fracture, nFractures, fractureSetSizeIncrement);
        if (nFractures == 0) initializeAttributeList();
        nFractures++;
        for (int i = 0; i < 4; i++)
        {
            StsPoint pt = new StsPoint(fracture.getCorner(i));
            addPoint(pt, 0, 0);
        }
    }

    public void addToModel()
    {
        StsProject project = currentModel.getProject();
        xOrigin = project.getXOrigin();
        yOrigin = project.getYOrigin();
        if(!project.addToProject(this, StsProject.TD_DEPTH))
            return;
        currentModel.checkAddToCursor3d(this);
        super.addToModel();
    }


    public void setIsVisible(boolean value)
    {
        isVisible = value;
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean getIsVisible()
    {
        return isVisible;
    }

    public void setStsColor(StsColor color)
    {
        stsColor = color;
        dbFieldChanged("stsColor", stsColor);
        colorscaleChanged = true;
        currentModel.viewObjectRepaint(this, this);
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }

    /**
     * Handle colorscale changes
     *
     * @param e - event
     */
    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            colorscaleChanged = true;
            currentModel.viewObjectChangedAndRepaint(this, this);
        }
    }

    /**
     * Set the colorscale to be applied to the colorby curve
     *
     * @param colorscale - selected colorscale
     */
    public void setColorscale(StsColorscale colorscale)
    { this.colorscale = colorscale; }

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            displayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsMVFractureSet.class, "isVisible", "Enable"),
                    //new StsBooleanFieldBean(StsMVFractureSet.class, "drawEdges", "Edges"),
                    totalFracsBean = new StsIntFieldBean(StsMVFractureSet.class, "fracCount", "Number Fractures:"),
                    fracsDisplayedBean = new StsIntFieldBean(StsMVFractureSet.class, "fracsDisplayed", "Enable Fractures:"),
                    new StsBooleanFieldBean(StsMVFractureSet.class, "displaySectionEdges", "Display Edges:"),
                    new StsColorComboBoxFieldBean(StsMVFractureSet.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors()),
                    decimationBean = new StsIntFieldBean(this, "decimationPercent", 0, 99, "% Decimation::", true),
                    verticalScaleBean = new StsFloatFieldBean(StsMVFractureSet.class, "verticalScale", true, "Vertical Scale",true),
                    colorByBean = new StsComboBoxFieldBean(StsMVFractureSet.class, "colorBy", "Color By:", "attributeList"),
                    colorMinBean = new StsFloatFieldBean(StsMVFractureSet.class, "colorMinimum", true, "Color Min"),
                    colorMaxBean = new StsFloatFieldBean(StsMVFractureSet.class, "colorMaximum", true, "Color Max"),
                    colorscaleBean = new StsEditableColorscaleFieldBean(StsMVFractureSet.class, "colorscale")
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
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet").selected(this);
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
}
