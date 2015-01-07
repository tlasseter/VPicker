package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 5:51:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFractureSetClass extends StsClass implements StsClassCursor3dTextureDisplayable, StsClassFractureDisplayable
{
    StsSpectrum spectrum = currentModel.getSpectrum("Basic");
    int nextColorIndex = 0;
    int newNameIndex = 0;

    private float defaultStimulatedRadius = 500.0f;
    private float gridSize = 0.1f;
    private float contourInterval = 0.0f;
    private byte displayType = DISPLAY_NONE;
    protected boolean displayOnSubVolumes = true;

    transient private StsRotatedGridBoundingBox fractureBoundingBox = null;
    transient private byte zDomain = StsParameters.NONE;
    transient public float rowGridSize;
    transient public float colGridSize;
    transient StsColorList distanceColorList;
    transient String displayTypeString = displayTypeStrings[displayType];

    static public final byte DISPLAY_NONE = 0;
    static public final byte DISPLAY_RADIUS = 1; // displays volume stimulated out to a radius away from fractures
    static public final byte DISPLAY_VOLUME = 2; // displays contour lines of distance away from fractures
    static final String[] displayTypeStrings = new String[] { "None", "Radius", "Volume" };
    static final int nBasicColors = StsFractureCursorSection.nBasicColors;
    public StsFractureSetClass()
    {
        userName = "Sets of Fractures";
    }

    public StsColor getNextColor()
    {
        return new StsColor(spectrum.getColor(nextColorIndex++));
    }

    public String getNextName()
    {
        return "fractureSet-" + newNameIndex++;
    }

    public StsCursor3dTexture constructDisplayableSection(StsModel model, StsCursor3d cursor3d, int dir)
    {
        return new StsFractureCursorSection(model, cursor3d, dir);
    }

    public boolean drawLast() { return true; }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
            {
                new StsComboBoxFieldBean(this, "displayTypeString", "Stimulation Display:", displayTypeStrings),
                new StsBooleanFieldBean(this, "displayOnSubVolumes", "Filter by SubVolumes")
            };
    }

    public void initializePropertyFields()
    {
        propertyFields = new StsFieldBean[]
            {
                new StsFloatFieldBean(this, "defaultStimulatedRadius", 0.0f, 10000.0f, "Default stimulated radius"),
                new StsFloatFieldBean(this, "gridSize", 0.0f, 1000.0f, "Display grid size"),
                new StsFloatFieldBean(this, "contourInterval", 0.0f, 1000.0f, "Display contour interval")
            };
    }

    public float getDefaultStimulatedRadius()
    {
        return defaultStimulatedRadius;
    }

    public void setDefaultStimulatedRadius(float defaultStimulatedRadius)
    {
        this.defaultStimulatedRadius = defaultStimulatedRadius;
        contourInterval = defaultStimulatedRadius/2;
    }

    public float getGridSize()
    {
        return gridSize;
    }

    public void setGridSize(float gridSize)
    {
        this.gridSize = gridSize;
    }

    public StsRotatedGridBoundingBox getFractureBoundingBox()
    {
        try
        {
            byte zDomainProject = currentModel.getProjectZDomain();
            if(fractureBoundingBox != null && zDomain == zDomainProject) return fractureBoundingBox;
            fractureBoundingBox = currentModel.getProject().getZDomainRotatedBoundingBox(zDomainProject);
            zDomain = zDomainProject;

            // define fractureBoundingBox xInc and yInc so that they are integral values of current xInc or yInc or vice versa
            float yInc = fractureBoundingBox.yInc;
            float ySize = fractureBoundingBox.getYSize();
            rowGridSize = computeIntegralGridSize(yInc, gridSize, ySize);
            int nRows = StsMath.floor(ySize/rowGridSize) + 1;

            float xInc = fractureBoundingBox.xInc;
            float xSize = fractureBoundingBox.getXSize();
            colGridSize = computeIntegralGridSize(xInc, gridSize, xSize);
            int nCols = Math.round(xSize/colGridSize) + 1;

            fractureBoundingBox.nRows = nRows;
            fractureBoundingBox.nCols = nCols;
            fractureBoundingBox.xInc = colGridSize;
            fractureBoundingBox.yInc = rowGridSize;
            float xMin = fractureBoundingBox.xMin;
            fractureBoundingBox.xMax = xMin + colGridSize*(nCols - 1);
            float yMin = fractureBoundingBox.yMin;
            fractureBoundingBox.yMax = yMin + rowGridSize*(nRows - 1);
            return fractureBoundingBox;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getFractureBoundingBox", e);
            return null;
        }
    }

    public float getRowGridSize()
    {
        checkFractureBoundingBox();
        return rowGridSize;
    }

    private boolean checkFractureBoundingBox()
    {
        return getFractureBoundingBox() != null;
    }

    public float getColGridSize()
    {
        checkFractureBoundingBox();
        return colGridSize;
    }

    static final int maxGridSize = 1000;

    /** compute a gridSize based on the gridSize provided such that it is an integral fraction of the current gridIncrement
     *  or an integral multiplier of the current gridIncrement subject to the limit of a maxGridSize
     * @param inc the current grid increment
     * @param gridSize desired grid cell size
     * @param gridLength length of the current grid in this direction
     * @return resulting grid cell size
     */
    private float computeIntegralGridSize(float inc, float gridSize, float gridLength)
    {
        int nGridCells = Math.round(gridLength/gridSize);
        if(nGridCells > maxGridSize)
        {
            float approxInc = Math.round(gridLength/maxGridSize);
            if(inc < approxInc)
            {
                int nIncsPerGrid = Math.round(approxInc/inc);
                return nIncsPerGrid*inc;
            }
            int nGridsPerInc = Math.round(inc/approxInc);
            return inc/nGridsPerInc;
        }
        float nGridsPerIncF = inc/gridSize;
        if(nGridsPerIncF > 1.0)
        {
            int nGridsPerInc = Math.round(nGridsPerIncF);
            return inc/nGridsPerInc;
        }
        else
        {
            int nIncsPerGrid = Math.round(1.0f/nGridsPerIncF);
            return inc*nIncsPerGrid;
        }
    }

    public boolean getDisplayOnSubVolumes() { return displayOnSubVolumes; }

    public void setDisplayOnSubVolumes(boolean b)
    {
        if (this.displayOnSubVolumes == b) return;
        this.displayOnSubVolumes = b;
//        setDisplayField("displayOnSubVolumes", displayOnSubVolumes);
        currentModel.subVolumeChanged();
        currentModel.win3dDisplayAll();
    }

    public float getBoundingBoxMaxSize()
    {
        //return 20*Math.max(rowGridSize, colGridSize);
        return Math.max(fractureBoundingBox.getXSize(), fractureBoundingBox.getYSize());
    }

    public boolean setColorList(GL gl, int shader, boolean displayOnOtherCursorSections)
	{
        float[][] arrayRGBA = computeRGBAArray(displayOnOtherCursorSections);
        return setColorList(gl, shader, arrayRGBA);
	}

    private boolean setColorList(GL gl, int shader, float[][] arrayRGBA)
     {
         if(arrayRGBA == null)return false;

         if(shader != StsJOGLShader.NONE)
         {
             StsJOGLShader.loadEnableARBColormap(gl, arrayRGBA, shader);
         }
         else
         {
             setColorList(gl, arrayRGBA);
         }
         arrayRGBA = null;
         return true;
     }

     /** If the number of fractureSets is n, there are n+3 colors.
      *  The first two colors are clear and black which make a screendoor transparency background if we are not displaying
      *  on another cursorSection.  If we are, the second color is also clear.
      *  The third color is the contour line color.
      *  The remaining colors are the colors of the fractureSets.
      */
     private float[][] computeRGBAArray(boolean displayOnOtherCursorSections)
     {
         int nFractureSets = getSize();
         int maxIndex = getLast().getIndex();
         int nColors = maxIndex + 3;
         nColors = StsMath.nextBaseTwoInt(nColors);
         float[][] arrayRGBA = new float[4][nColors];
         if(displayOnOtherCursorSections)
         {
            setRGBA(StsColor.BLACK, 1, arrayRGBA, 1.0f);
         }
         else
         {
//            setRGBA(StsColor.BLACK, 0, arrayRGBA, 0.5f); // make plane background semitransparent
            setRGBA(StsColor.WHITE, 1, arrayRGBA, 1.0f);
         }
         for(int n = 0; n < nFractureSets; n++)
         {
             StsFractureSet fractureSet = (StsFractureSet)getElement(n);
             StsColor color = fractureSet.getStsColor();
             int index = fractureSet.getIndex();
             setRGBA(color, index+nBasicColors, arrayRGBA, 0.5f);
         }
         return arrayRGBA;
     }

     private void setRGBA(StsColor color, int index, float[][] arrayRGBA, float alpha)
     {
        float[] rgba = color.getRGBA();
        for (int i = 0; i < 3; i++)
            arrayRGBA[i][index] = rgba[i];
        arrayRGBA[3][index] = alpha;
     }

     private void setColorList(GL gl, float[][] arrayRGBA)
     {
         int nColors = arrayRGBA[0].length;
         gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
         gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
         gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
         gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
         gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
//        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
         arrayRGBA = null;
     }

    public Iterator getFractureIterator(boolean visibleOnly)
    {
        return new FractureIterator(visibleOnly);
    }

    public byte getDisplayType() { return displayType; }
    public void setDisplayType(byte displayType)
    {
        this.displayType = displayType;
        displayTypeString = displayTypeStrings[displayType];
    }
    public String getDisplayTypeString() { return displayTypeStrings[displayType]; }
    public void setDisplayTypeString(String displayTypeString)
    {
        for(int n = 0; n < 3; n++)
        {
            if(displayTypeString == displayTypeStrings[n])
            {
                this.displayTypeString = displayTypeString;
                byte newDisplayType = (byte)n;
                if(this.displayType == newDisplayType) return;
                this.displayType = newDisplayType;
                currentModel.viewObjectChangedAndRepaint(this, getFirst());
            }
        }
    }

    public void setDisplayTypeNone() { setDisplayType(DISPLAY_NONE); }
    public void setDisplayTypeRadius() { setDisplayType(DISPLAY_RADIUS); }
    public void setDisplayTypeVolume() { setDisplayType(DISPLAY_VOLUME); }

    public boolean isDisplayingVolume() { return displayType != DISPLAY_NONE; }
    public boolean isDisplayingRadius() { return displayType == DISPLAY_RADIUS; }
    public boolean displayStimulatedVolume() { return displayType != DISPLAY_NONE; }

    public float getContourInterval()
    {
        return contourInterval;
    }

    public void setContourInterval(float contourInterval)
    {
        this.contourInterval = contourInterval;
    }

    /** This iterates over all fractures in all fractureSets */
    class FractureIterator implements Iterator
    {
        boolean visibleOnly;
        int nFractureSet = 0;
        StsFractureSet fractureSet = null;
        Iterator fractureSetIterator;
        Iterator fractureIterator;

        FractureIterator(boolean visibleOnly)
        {
            this.visibleOnly = visibleOnly;
            fractureSetIterator = getObjectIterator(visibleOnly);

        }

        public boolean hasNext()
        {
            if(fractureIterator == null || !fractureIterator.hasNext())
            {
                if(fractureSetIterator.hasNext())
                {
                    fractureSet = (StsFractureSet)fractureSetIterator.next();
                    fractureIterator = fractureSet.getFractureList().getIterator();
                }
                else
                    return false;
            }
            return fractureIterator.hasNext();
        }

        public Object next()
        {
            return fractureIterator.next();
        }

        public void remove() { }        
    }

	public ArrayList<StsFractureDisplayable> getDisplayableFractures()
	{
		ArrayList<StsFractureDisplayable> displayables = new ArrayList<StsFractureDisplayable>();
		int nFractureSets = getSize();
		for(int n = 0; n < nFractureSets; n++)
		{
			StsFractureSet fractureSet = (StsFractureSet)getElement(n);
			StsObjectRefList fractureList = fractureSet.getFractureList();
		    int nFractures = fractureList.getSize();
			for(int i = 0; i < nFractures; i++)
				displayables.add((StsFractureDisplayable)fractureList.getElement(n));
		}
		return displayables;
	}
}