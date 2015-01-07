package com.Sts.DBTypes;


import com.Sts.Actions.Wizards.PostStack2dLoad.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 22, 2007
 * Time: 4:06:59 PM
 * To change this template use File | Settings | File Templates.
 */

public class StsSeismicLineSet extends StsSeismic implements StsTreeObjectI
{
	public StsSeismicLine2d[] lines;
    public StsSeismicLineSetClass lineSetClass;
    float zShift = 0.0f;

    static public StsFieldBean[] displayFields2d = null;
	static private StsFieldBean[] propertyFields2d = null;
    static protected StsObjectPanel objectPanel = null;

    static StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsSeismicLineSet.class, "colorscale");

    static protected StsComboBoxFieldBean lineListBean = new StsComboBoxFieldBean(StsSeismicLineSet.class, "line", "Current Line:");
    static protected StsFloatFieldBean zShiftBean = new StsFloatFieldBean(StsSeismicLineSet.class, "zShift", -1000, 1000, "Z Shift:", true);
    transient StsSeismicLine2d currentLine = null;

    public String getGroupname()
    {
        return groupNone;
    }

    public StsSeismicLineSet()
	{
	}

	public StsSeismicLineSet(boolean persistent)
	{
		super(persistent);
        lineSetClass = (StsSeismicLineSetClass) getCreateStsClass();
    }

    public boolean initialize(StsModel model)
    {
        lineSetClass = (StsSeismicLineSetClass) getCreateStsClass();
        initializeColorscale();
        lineListBean.setListItems(lines);
        setLine(lines[0]);
        return true;
    }

    public StsSeismicLineSetClass getSeismicLineSetClass()
	{
		return (StsSeismicLineSetClass)currentModel.getCreateStsClass(getClass());
	}

    public void addToModel()
    {
        super.addToModel();
        if(lines == null) return;
        for(int n = 0; n < lines.length; n++)
        {
            lines[n].addToModel();
        }
        currentModel.getProject().checkAddUnrotatedClass(StsSeismicLineSet.class);
    }

    public void setLinesXYs()
    {
        for(int n = 0; n < lines.length; n++)
            lines[n].setLineXYsFromCDPs();
    }


	public void setDataHistogram()
	{
		if (dataHist != null && colorscaleBean != null)
			colorscaleBean.setHistogram(dataHist);
	}

    public void display(StsGLPanel3d glPanel)
	{
        if(!isVisible)
            return;
		GL gl = glPanel.getGL();
		if (lines != null && lines.length > 0)
		{
//			checkMemoryStatus("before display 2d lines");
			for (int i = 0; i < lines.length; i++)
            {
				lines[i].display(glPanel);
            }
//			checkMemoryStatus("after display 2d lines");
			gl.glEnable(GL.GL_LIGHTING);
		}
	}
    public float getZShift()
    {
        return currentLine.getZShift();
    }
    public void setZShift(float zshift)
    {
        currentLine.setZShift(zshift);
        currentModel.viewObjectRepaint(this, this);
    }

    public void actionPerformed(ActionEvent e)
    {
        if (!(e.getSource() instanceof StsColorscale))
            return;
        seismicColorList.setColorListChanged(true);
        for(int i=0; i<lines.length; i++)
        {
            lines[i].textureChanged(true);
        }
        currentModel.viewObjectRepaint(this, this);
    }

    public boolean textureChanged()
    {
       for(int i=0; i<lines.length; i++)
            lines[i].textureChanged();
       return true;
    }

    //TODO implement this
    public boolean getTraceValues(int row, int col, int sliceMin, int sliceMax, int dir, boolean useByteCubes, float[] floatData)
    {
        return false;
    }

    public int getIntValue(int row, int col, int slice)
    {
        return 0;
    }

    public StsFieldBean[] getDisplayFields()
    {
        try
        {
            if (displayFields2d == null)
            {
                colorscaleBean = new StsEditableColorscaleFieldBean(StsSeismicLineSet.class, "colorscale");
                zShiftBean.fixStep(zInc);
                displayFields2d = new StsFieldBean[]
                                  {
                                          new StsBooleanFieldBean(StsSeismicLineSet.class, "isVisible", "Enable"),
                                          lineListBean, zShiftBean, colorscaleBean
                                  };
            }
            colorscaleBean.setValueObject(colorscale);
            colorscaleBean.setHistogram(dataHist);
            return displayFields2d;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSeismicLineSet.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void setLine(StsSeismicLine2d line)
    {
        currentLine = line;
        zShiftBean.setValue(currentLine.getZShift());
    }
    public StsSeismicLine2d getLine() { return currentLine; }

    public StsFieldBean[] getPropertyFields()
	{
		try
		{
			if (propertyFields2d == null)
			{
				propertyFields2d = new StsFieldBean[]
					{
					new StsStringFieldBean(StsSeismicLineSet.class, "name", true, "Name"),
					new StsStringFieldBean(StsSeismicLineSet.class, "zDomainString", false, "Z Domain"),
					new StsDoubleFieldBean(StsSeismicLineSet.class, "xOrigin", false, "X Origin"),
					new StsDoubleFieldBean(StsSeismicLineSet.class, "yOrigin", false, "Y Origin")
				};
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeismicLineSet.getPropertyFields() failed.", e, StsException.WARNING);
			return null;
		}
		return propertyFields2d;
	}

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        lineSetClass.selected(this);
        currentModel.getGlPanel3d().checkAddView(StsView3d.class);
//        currentModel.glPanel3d.cursor3d.setCurrentSeismicVolume(this);
        currentModel.win3dDisplayAll();
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public void initializeHistogram(float excludeValue)
    {
        int dataCnt[] = new int[255];
        int ttlSamples = 0;
        for(int i=0; i<lines.length; i++)
        {
            float[] data = lines[i].readFloatData();
            if(data == null) continue;
            ttlSamples =+ data.length;
            for(int j=0; j<data.length; j++)
            {
                if(data[j] != excludeValue)
                    StsToolkit.accumulateHistogram(dataCnt, dataMin, dataMax, data[j]);
            }
        }
        dataHist = StsToolkit.calculateHistogram(dataCnt,ttlSamples);
    }

    public boolean add2dLines(StsFile files[], StsProgressPanel panel, StsLine2dLoad step)
	{
		if (files == null) return false;
        int nFiles = files.length;
        int nExistingFiles = 0;
        for (int n = nFiles-1; n >= 0; n--)
        {
            if (!files[n].exists())
            {
                panel.setDescription(": File not found " + files[n].getPathname());
                StsMath.arrayDeleteElement(files, n);
            }
            else
                nExistingFiles++;
        }
        nFiles = nExistingFiles;
        lines = new StsSeismicLine2d[nFiles];
        double progressStep = 1.0 / (nFiles + 1);
        int nn = 0;
        byte lineSetZDomain = StsProject.TD_NONE;
        panel.initialize(nFiles);
        for (int n = 0; n < nFiles; n++)
        {
            StsSeismicLine2d line = StsSeismicLine2d.constructor(files[n], currentModel, this);
            if(line == null) continue;
            if (nn == 0)
            {
                setName(line.stemname);
                stsDirectory = line.stsDirectory;
                lineSetZDomain = line.getZDomain();
                setZDomain(lineSetZDomain);
            }
            else
            {
                if(lineSetZDomain != line.getZDomain())
                {
                    new StsMessage(currentModel.win3d, StsMessage.WARNING, "Loading only lines in " + StsParameters.TD_ALL_STRINGS[lineSetZDomain] + ".\n" +
                            " So will not load line " + line.stemname + " which is " + StsParameters.TD_ALL_STRINGS[line.getZDomain()]);
                    continue;
                }
            }
            lines[nn++] = line;
            dataMin = Math.min(dataMin, line.dataMin);
            dataMax = Math.max(dataMax, line.dataMax);
            addUnrotatedBoundingBox(line);
            panel.setValue(n+1);
            panel.setDescription("Loaded line : " + (n+1) + " of " + nFiles);
            panel.appendLine("Completed loading line index:" + line.getIndex() + " from file " + line.getName());
        }

        StsProject project = currentModel.getProject();
        if(!project.addToProject(this, false))
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to add lines to volume: " + files[0].getPathname());
            return false;
        }

        lines = (StsSeismicLine2d[])StsMath.trimArray(lines, nn);
        setLinesXYs();
        setLinesIndexRanges(lines);

        currentModel.setCurrentObject(this);
        lineListBean.setListItems(lines);
        setLine(lines[0]);
        return true;
    }
}
