package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;

public class StsPlatform extends StsMainTimeObject implements StsSerializable, StsTreeObjectI
{
    int nSlots = 16;
    int numRows = 4;
    int numCols = 4;
    float rowSpacing = 8.0f;
    float colSpacing = 8.0f;
//    boolean isVisible = true;
    String[] slotNames = null;
    String[] wellNames = null;
//    int[] wellIndices = null;
    double xOrigin = nullValue;
    double yOrigin = nullValue;
    float zKB;
//    StsPoint origin = new StsPoint(0.0f, 0.0f, 0.0f);
    float rotationAngle = 0.0f;
    byte slotType = RECTANGULAR;
    byte platformType = UNKNOWN;
    /** The color of this platform described in red-green-blue-alpha values from 0 to 1. */
    protected StsColor stsColor = new StsColor(StsColor.BLUE);

    transient JPanel slotGraph = null;
    transient int currentSlot = -1;
    static StsObjectPanel objectPanel = null;

    transient public static final byte RECTANGULAR = 0;
    transient public static final byte CIRCULAR = 1;
    static public final String[] SLOT_LAYOUT_STRINGS = new String[] { "Rectangular", "Circular" };
    static public final byte[] SLOT_LAYOUT = new byte[] { RECTANGULAR, CIRCULAR };

    transient public static final byte UNKNOWN = 0;
    transient public static final byte FIXED = 1;
    transient public static final byte JACKUP = 2;
    transient public static final byte SEMISUB = 3;
    transient public static final byte TENSIONLEG = 4;
    transient public static final byte DRILLSHIP = 5;
    static public final String[] PLATFORM_TYPE_STRINGS = new String[] { "Unknown", "Fixed" , "Jackup", "Semi-Submersible", "Tension Leg", "Drill Ship" };
    static public final byte[] PLATFORM_TYPE = new byte[] { UNKNOWN, FIXED, JACKUP, SEMISUB, TENSIONLEG, DRILLSHIP };

    static StsDateFieldBean bornField = null;
    static StsDateFieldBean deathField = null;

    static final double nullValue = StsParameters.nullValue;

	static final long serialVersionUID = 1L;

    static public StsFieldBean[] displayFields = null;

    public StsPlatform(boolean persistent)
    {
        super(persistent);
        setName("Platform");
        configurePlatform();
    }

    public StsPlatform()
    {
    }

    public StsPlatform(String name, int nRows, int nCols, boolean persistent)
    {
        super(persistent);
        setName(name);
        numRows = nRows;
        numCols = nCols;
        configurePlatform();
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            super.isVisible = isVisible;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPlatform.classInitialize() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public void configurePlatform()
    {
        if(slotNames != null)
        {
            if((numCols * numRows) > nSlots)
            {
                slotNames = (String[]) StsMath.arrayGrow(slotNames,(numCols * numRows) - nSlots);
                wellNames = (String[]) StsMath.arrayGrow(wellNames,(numCols * numRows) - nSlots);
                for(int i=nSlots-1; i<numCols * numRows; i++)
                {
                    slotNames[i] = String.valueOf(i + 1);
                    wellNames[i] = null;
                }
            }
            else
            {
                slotNames = (String[])StsMath.trimArray(slotNames, numCols * numRows);
                wellNames = (String[])StsMath.trimArray(wellNames, numCols * numRows);
            }
        }
        else
        {
            slotNames = new String[(numCols * numRows)];
            wellNames = new String[(numCols * numRows)];
            for(int i=0; i<(numCols * numRows); i++)
            {
                slotNames[i] = String.valueOf(i + 1);
                wellNames[i] = null;
            }
        }
        nSlots = numCols * numRows;
        StsPlatformClass platformClass = (StsPlatformClass)currentModel.getStsClass(StsPlatform.class);
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
		if (platformClass != null && spectrum != null)
		{
			stsColor = platformClass.getNextColor(spectrum);
		}
    }

    public String getWellAtSlot(String slotName)
    {
        int slotIdx = getIndexFromSlotName(slotName);
        if(slotIdx != -1)
        {
            StsWell well = (StsWell)currentModel.getObjectWithName(StsWell.class, wellNames[slotIdx]);
            if(well != null) return wellNames[slotIdx];
            StsWellPlan plannedWellSet = (StsWellPlan)currentModel.getObjectWithName(StsWell.class, wellNames[slotIdx]);
            if(plannedWellSet != null) return wellNames[slotIdx];

        }
        return null;
    }

    public float getWellDirectionAngle(String slotName)
    {
        int slotIdx = getIndexFromSlotName(slotName);
        if (slotIdx != -1)
        {
            StsWell well = (StsWell) currentModel.getObjectWithName(StsWell.class, wellNames[slotIdx]);
            if (well != null)
                return well.getWellDirectionAngle();
            StsWellPlan plannedWellSet = (StsWellPlan) currentModel.getObjectWithName(StsWellPlan.class, wellNames[slotIdx]);
            if (plannedWellSet != null)
                return plannedWellSet.getWellDirectionAngle();
        }
        return StsParameters.nullValue;
    }
    public void clearSlotAssignment()
    {
        for(int i=0; i<nSlots; i++)
            wellNames[i] = null;
        drawConfiguration();
    }

    public void clearSlotAssignment(int slotIdx)
   {
       wellNames[slotIdx] = null;
       drawConfiguration();
   }

   public void clearSlotAssignment(String wellname)
  {
      int slotIdx = this.getIndexFromWellName(wellname);
      if(slotIdx != -1)
          wellNames[slotIdx] = null;
      drawConfiguration();
  }

    public boolean addWellWithSlotName(String slotName, String wellname)
    {
        int slotIdx = getIndexFromSlotName(slotName);
        return addWellAtSlot(slotIdx, wellname);
    }

    public boolean addWellAtSlot(int slotIdx, String wellname)
    {
        StsPlatformClass pc = (StsPlatformClass)currentModel.getStsClass(StsPlatform.class);
        if(slotIdx != -1)
        {
            if((wellNames[slotIdx] != null) && (!wellname.equals(wellNames[slotIdx])))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING,
                               "Well " + wellNames[slotIdx] + " already assigned to slot " +  slotNames[slotIdx]);
                return false;
            }
            else if((pc.getWellPlatform(wellname) != null) && (pc.getWellPlatform(wellname) != this))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING,
                               "Well " + wellname + " already assigned to slot " +  pc.getWellPlatform(wellname).getSlotName(wellname)
                               + " on platform " + pc.getWellPlatform(wellname).getName());
                return false;
            }
            // Assign well to slot
            wellNames[slotIdx] = wellname;
            this.dbFieldChanged("wellNames", wellNames);
            return true;
        }
        return false;
    }

    public boolean deleteWellfromPlatform(String wellname)
    {
        int slotIdx = getIndexFromWellName(wellname);
        if(slotIdx > -1)
        {
            wellNames[slotIdx] = null;
            if(getNWells() == 0)
                delete();
            else
                dbFieldChanged("wellNames", wellNames);
            return true;
        }
        return false;
    }

    public int getIndexFromWellName(String wellName)
    {
		if(wellName == null || wellNames == null) return -1;
        for(int i=0; i<wellNames.length; i++)
        {
            if (wellName.equals(wellNames[i]))
            {
                return i;
            }
        }
        return -1;
    }

    public int getNWells()
    {
        int nWells = 0;
        for(int i = 0; i < wellNames.length; i++)
            if(wellNames[i] != null) nWells++;
        return nWells;
    }

    public int getIndexFromSlotName(String name)
    {
        for(int i=0; i<slotNames.length; i++)
        {
            if(name.equals(slotNames[i]))
            {
                return i;
            }
        }
        return -1;
    }

    public boolean deleteWellAtSlot(String slotName)
    {
        int slotIdx = getIndexFromSlotName(slotName);
        return deleteWellAtSlot(slotIdx);
    }

    public boolean deleteWellAtSlot(int slotIdx)
    {
        if(slotIdx != -1)
        {
            wellNames[slotIdx] = null;
            return true;
        }
        return false;
    }

    public boolean isWellInSlot(int idx)
    {
        if(wellNames[idx] != null)
            return true;
        else
            return false;
    }

    public float getRotationAngle() { return rotationAngle; }
    public void setRotationAngle(float angle)
    {
        rotationAngle = angle;
        drawConfiguration();
    }
    public double getXOrigin() { return xOrigin; }
    public void setXOrigin(double x)
    {
        xOrigin = x;
//        origin.setX(x);
//        drawConfiguration();
    }
    public double getYOrigin() { return yOrigin; }
    public void setYOrigin(double y)
    {
        yOrigin = y;
//        origin.setY(y);
//        drawConfiguration();
    }
    public float getZKB() { return zKB; }
    public void setZKB(float z)
    {
        zKB = z;
//        origin.setZ(z);
//        drawConfiguration();
    }
    public int getNCols() { return numCols; }
    public void setNCols(int nCols)
    {
        numCols = nCols;
        configurePlatform();
        drawConfiguration();
    }
    public float getColSpacing() { return colSpacing; }
    public void setColSpacing(float spacing)
    {
        colSpacing = spacing;
        drawConfiguration();
    }

    public int getNRows() { return numRows; }
    public void setNRows(int nRows)
    {
        numRows = nRows;
        configurePlatform();
        drawConfiguration();
    }
    public float getRowSpacing() { return rowSpacing; }
    public void setRowSpacing(float spacing)
    {
        rowSpacing = spacing;
        drawConfiguration();
    }
    public int getNumSlots() { return nSlots; }
    public String getSlotName(String wellname) { return slotNames[getIndexFromWellName(wellname)]; }
//    public StsPoint getOrigin() { return origin; }
//    public void setOrigin(StsPoint point) { origin = point; }
    public String[] getSlotNames() { return slotNames; }
    public void setSlotNames(String[] names)
    {
        if(names.length == slotNames.length)
            slotNames = names;
        else
            new StsMessage(currentModel.win3d, StsMessage.WARNING, "Number of names must equal number of slots(" + nSlots + ")");
    }
    public String getSlotName(int idx) { return slotNames[idx]; }
    public void setSlotName(int idx, String name) { slotNames[idx] = name; }
    public Object[] getWellNames() { return wellNames; }
    public boolean isWellOnPlatform(String wellname)
    {
        if(getIndexFromWellName(wellname) != -1)
            return true;
        return false;
    }

    public void setName(String name)
    {
        super.setName(name);
        drawConfiguration();
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setStsColor(StsColor color)
    {
        if(stsColor.equals(color)) return;
        stsColor = color;
        currentModel.win3dDisplayAll();
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;
		bornField = new StsDateFieldBean(StsPlatform.class, "bornDate", "Born Date:");
		deathField = new StsDateFieldBean(StsPlatform.class, "deathDate", "Death Date:");
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(StsPlatform.class, "isVisible", "Enable:"),
			bornField,
			deathField,
            new StsStringFieldBean(StsPlatform.class, "name", true, "Name:"),
            new StsStringFieldBean(StsPlatform.class, "platformTypeByString", false, "Type:"),
            new StsFloatFieldBean(StsPlatform.class, "xOrigin", false, "X Origin:"),
            new StsFloatFieldBean(StsPlatform.class, "yOrigin", false, "Y Origin:"),
            new StsFloatFieldBean(StsPlatform.class, "zKB", false, "KB Depth:"),
            new StsFloatFieldBean(StsPlatform.class, "rotationAngle", false, "Rotation Angle:"),
            new StsStringFieldBean(StsPlatform.class, "slotLayoutByString", false, "Layout:"),
            new StsIntFieldBean(StsPlatform.class, "numSlots", false, "Number of Slots:"),
            new StsIntFieldBean(StsPlatform.class, "numRows", false, "Number of Rows:"),
            new StsFloatFieldBean(StsPlatform.class, "rowSpacing", false, "Row Spacing:"),
            new StsIntFieldBean(StsPlatform.class, "numCols", false, "Number of Columns:"),
            new StsFloatFieldBean(StsPlatform.class, "colSpacing", false, "Column Spacing:"),
            new StsColorComboBoxFieldBean(StsPlatform.class, "stsColor", "Color")
       };
       return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
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

    public boolean anyDependencies()
    {
        return false;
    }

    public void treeObjectSelected()
    {
        ;
    }

    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        isVisible = b;
        currentModel.win3dDisplayAll();
    }
    public boolean getIsVisible() { return isVisible; }

    public void updatePanel()
    {
        ;
    }
    public String getPlatformTypeByString()
    {
        return PLATFORM_TYPE_STRINGS[platformType];
    }
    public void setPlatformTypeByString(String typeString)
    {
        for(int i=0; i<PLATFORM_TYPE_STRINGS.length; i++)
        {
            if (typeString.equals(PLATFORM_TYPE_STRINGS[i]))
            {
                platformType = PLATFORM_TYPE[i];
                break;
            }
        }
        return;
    }

    public String getSlotLayoutByString()
    {
        return SLOT_LAYOUT_STRINGS[slotType];
    }
    public void setSlotLayoutByString(String typeString)
    {
        for(int i=0; i<SLOT_LAYOUT_STRINGS.length; i++)
        {
            if (typeString.equals(SLOT_LAYOUT_STRINGS[i]))
            {
                slotType = SLOT_LAYOUT[i];
                drawConfiguration();
                break;
            }
        }
        return;
    }

    public void setCanvas(JPanel graph)
    {
        slotGraph = graph;
    }

    public void drawConfiguration()
    {
        if(slotGraph == null)
            return;
        slotGraph.repaint();
    }

    public boolean checkValidity()
    {
        if((colSpacing < 1.0) || (rowSpacing < 1.0))
            return false;
        if((numCols < 1) || (numRows < 1))
            return false;
        return true;
    }

    public void setCurrentSlotIndex(int idx)
    {
        currentSlot = idx;
    }

    public int getCurrentSlotIndex() { return currentSlot; }

    public boolean addWellToCurrentSlot(String wellname)
    {
        if(currentSlot == -1) return false;
        return addWellAtSlot(currentSlot, wellname);
    }

    public void deleteWellAtCurrentSlot()
    {
        if(currentSlot == -1) return;
        this.deleteWellAtSlot(currentSlot);
    }

    public double[] getSlotXY()
    {
        double rowCenter = (float)(numRows-1)/2;
        double colCenter = (float)(numCols-1)/2;
        double row = currentSlot/numCols;
        double col = currentSlot - row*numCols;
        double rowOffset = row - rowCenter;
        double colOffset = col - colCenter;
        double y = rowOffset*rowSpacing;
        double x = colOffset*colSpacing;
        double angleRad = rotationAngle*StsMath.RADperDEG;
        double cosA = Math.cos(angleRad);
        double sinA = Math.sin(angleRad);
        double rotatedY = y*cosA - x*sinA;
        double rotatedX = x*cosA + y*sinA;
        return new double[] { xOrigin + rotatedX, yOrigin + rotatedY };
    }

    public void display(StsGLPanel3d glPanel3d, float halfSize)
    {
        GL gl = glPanel3d.getGL();
        try
        {
            float[] xy = currentModel.getProject().getRelativeXY(xOrigin, yOrigin);
            float x = xy[0];
            float y = xy[1];
            float z = -zKB;
            if (!isDepth)
            {
                float estimatedT = currentModel.getProject().getTimeFromDepth(z);
               float[] xyzmt = new float[] { x, y, z, 0.0f, estimatedT };
                   StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
                   if (velocityModel != null)
                   z = (float) velocityModel.getT(xyzmt);
               else
                   z = estimatedT;
            }
            float templateXSize = colSpacing * numCols / 2;
            float templateYSize = rowSpacing * numRows / 2;
            float padXSize = Math.min(2 * templateXSize, halfSize);
            float padYSize = Math.min(2 * templateYSize, halfSize);
            drawPlatformRectangle(gl, x, y, z, 2 * padXSize, 2 * padYSize, stsColor);
            glPanel3d.setViewShift(gl, 1.0);
            drawPlatformRectangle(gl, x, y, z, templateXSize, templateYSize, StsColor.GREY);
        }
        catch(Exception e)
        {
            StsException.outputException("StsPlatform.display() failed.", e, StsException.WARNING);
        }
        finally
        {
            glPanel3d.resetViewShift(gl);
        }
    }

    private void drawPlatformRectangle(GL gl, float x, float y, float z, float halfXSize, float halfYSize, StsColor stsColor)
    {
        stsColor.setGLColor(gl);
        gl.glBegin(GL.GL_QUADS);
        if(this.rotationAngle == 0.0f)
        {
            gl.glVertex3f(x - halfXSize, y - halfYSize, z);
            gl.glVertex3f(x + halfXSize, y - halfYSize, z);
            gl.glVertex3f(x + halfXSize, y + halfYSize, z);
            gl.glVertex3f(x - halfXSize, y + halfYSize, z);
        }
        else
        {
            float cos = (float)StsMath.cosd(rotationAngle);
            float sin = (float)StsMath.sind(rotationAngle);
            gl.glVertex3f(x + halfXSize*(- cos - sin), y + halfYSize*(- cos + sin), z);
            gl.glVertex3f(x + halfXSize*(+ cos - sin), y + halfYSize*(- cos - sin), z);
            gl.glVertex3f(x + halfXSize*(+ cos + sin), y + halfYSize*(+ cos - sin), z);
            gl.glVertex3f(x + halfXSize*(- cos + sin), y + halfYSize*(+ cos + sin), z);
        }
        gl.glEnd();
        gl.glEnable(GL.GL_LIGHTING);
    }

    public boolean isOriginSet() { return xOrigin != nullValue && yOrigin != nullValue; }
    public void setBornDate(String born)
    {
        if(!StsDateFieldBean.validateDateInput(born))
        {
            bornField.setValue(StsDateFieldBean.convertToString(bornDate));
            return;
        }
        super.setBornDate(born);
    }
    public void setDeathDate(String death)
    {
        if(!StsDateFieldBean.validateDateInput(death))
        {
            deathField.setValue(StsDateFieldBean.convertToString(deathDate));
            return;
        }
        super.setDeathDate(death);
    }
}
