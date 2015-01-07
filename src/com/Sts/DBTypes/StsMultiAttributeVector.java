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
public class StsMultiAttributeVector extends StsMainTimeObject implements StsSelectable, StsTreeObjectI
{
    public static final byte LENGTH = 0;
    public static final byte AZIMUTH = 1;
    public static final byte COLOR = 2;

    public static final byte SOLID = 0;
    public static final byte DASHED = 1;

    public static final byte LINE = 0;
    public static final byte LINEENDS = 1;
    public static final byte DIAMOND = 2;
    public static final byte BOWTIE = 3;
    public static String[] typeStrings = { "LINE", "LINE & ENDPOINTS", "DIAMOND", "BOWTIE" };
    static StsObjectPanel objectPanel = null;

    static public StsFieldBean[] displayFields = null;
    static public StsComboBoxFieldBean cultureDisplayableList;
    static public StsFieldBean[] propertyFields = null;

    public StsCultureDisplayable cultureDisplayableObject = null;
    protected StsSeismicBoundingBox[] volumes = null;
    protected boolean normalizeAzimuth = true;
    protected float lengthThreshold = StsParameters.nullValue;
    protected boolean applyThreshold = false;
    protected boolean drawMidpoints = false;
    protected float dataMin = StsParameters.largeFloat;
    protected float dataMax = -StsParameters.largeFloat;
    protected int decimation = 10;
    protected int scale = 2;
    protected int width = 5;
    protected StsColorscale colorscale = null;
    protected String spectrumName = "Basic";
    protected byte symbolType = LINE;

    int stroke = SOLID;
    transient int numRows = 0, numCols = 0;
    transient boolean isPlanar = true;
    transient byte[][] planes = new byte[3][];
    transient float currentZ = 0.0f;
    transient StsCultureDisplayable currentCultureDisplayable = null;
    transient boolean usingDisplayLists = true;
    transient private int gridDisplayListNum = 0;

    transient boolean displayListChanged = true;
    
    static protected StsEditableColorscaleFieldBean colorscaleBean = new StsEditableColorscaleFieldBean(StsMultiAttributeVector.class, "colorscale");

    public StsMultiAttributeVector()
    {
    }

    public StsMultiAttributeVector(boolean persist)
    {
        setStsColor(StsColor.RED);
    }

    public StsMultiAttributeVector(String name, StsSeismicBoundingBox[] vols, boolean normalizeAz, float threshold)
    {
        super(false);
        setName(name);
        this.volumes = vols;
        cultureDisplayableObject = currentModel.getProject();
        dataMin = vols[COLOR].getDataMin();
        dataMax = vols[COLOR].getDataMax();
        initializeColorscale();
        colorscaleBean.setHistogram(((StsSeismicVolume)vols[COLOR]).dataHist);
        normalizeAzimuth = normalizeAz;
        lengthThreshold = threshold;
        if(lengthThreshold != StsParameters.nullValue)
            applyThreshold = true;
        initialize();
    }

    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
        	clearDisplayList();
            currentModel.viewObjectRepaint(this, this);
        }
        return;
    }

    static public StsMultiAttributeVectorClass getMultiAttributeVectorClass()
    {
        return (StsMultiAttributeVectorClass) currentModel.getCreateStsClass(StsMultiAttributeVector.class);
    }

    private void initializeColorscale()
    {
    	if(colorscale == null)
    	{
    		spectrumName = getMultiAttributeVectorClass().getDefaultSpectrumName();
    		colorscale = new StsColorscale("VectorSet", currentModel.getSpectrum(spectrumName),(float)dataMin, (float)dataMax);
    	}
    	colorscale.addActionListener(this);
    }

    public StsColorscale getColorscale() { return colorscale; }
    public void setColorscale(StsColorscale colorscale) { this.colorscale = colorscale; }
    public void initialize() { }
    public boolean initialize(StsModel model)
    {
        cultureDisplayableObject = model.getProject();
        initializeColorscale();
        colorscaleBean.setHistogram(((StsSeismicVolume)volumes[COLOR]).dataHist);
        if(lengthThreshold != StsParameters.nullValue)
            applyThreshold = true;
        return true;
    }

    public void setVolumes(StsSeismicVolume[] vols)
    {
        volumes = vols;
    }

    public boolean checkLoadPlanes(GL gl)
    {
        boolean validPlane = true;
        float z = 0.0f;

        if(currentCultureDisplayable == cultureDisplayableObject)
        {
        	if(isPlanar)
        	{
        		if(currentZ == currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR))
        			return true;
        	}
        	else
        	{
        		if(currentZ == ((StsSurface)cultureDisplayableObject).getOffset())
        			return true;
        	}
        }
        
        clearDisplayList(gl);
        
        isPlanar = cultureDisplayableObject.isPlanar();
        currentCultureDisplayable = cultureDisplayableObject;
        if(isPlanar)
        {   	
        	currentZ = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
        	for(int i = 0 ; i <volumes.length; i++)
        	{
        		if(volumes[i] == null)
        		{
        			StsMessageFiles.infoMessage("Unable to read one of the MVA volume planes.");
        			continue;
        		}
        		planes[i] = ((StsSeismicVolume)volumes[i]).readBytePlaneData(StsCursor3d.ZDIR, currentZ);
        		if(planes[i] == null)
        		{
        			validPlane = false;
        			StsMessageFiles.infoMessage("Unable to read one of the MVA volume planes.");
        			continue;
        		}
        	}
        	if(volumes[0] != null)
        	{
        		numRows = volumes[0].nRows;
        		numCols = volumes[0].nCols;
        	}
        }
        else // cultureDisplayableObject is a surface
        {
            StsSurface surface = (StsSurface)cultureDisplayableObject;
            for(int i = 0 ; i <volumes.length; i++)
            {
                if(volumes[i] == null)
                {
                    StsMessageFiles.infoMessage("Unable to read one of the MVA volume planes.");
                    continue;
                }
                currentZ = surface.getOffset();
                planes[i] = ((StsSeismicVolume)volumes[i]).getSurfaceTexture(surface).getTextureData();
                if(planes[i] == null)
                {
                    validPlane = false;
                    StsMessageFiles.infoMessage("Unable to read one of the MVA volume planes.");
                    continue;
                }
            }
            if(volumes[0] != null)
            {
                numRows = ((StsSurface)cultureDisplayableObject).nRows;
                numCols = ((StsSurface)cultureDisplayableObject).nCols;
            }
        }
        if(!validPlane)
            return false; 
        else
        	return true;
    }
    
    private void clearDisplayList()
    {
    	displayListChanged = true;
    }
    private void clearDisplayList(GL gl)
    {
        if (gridDisplayListNum > 0)
        {
            gl.glDeleteLists(gridDisplayListNum, 1);
            gridDisplayListNum = 0;
        }    	
    }
    /**
     * draw
     */
    public void display2d(StsGLPanel glPanel)
    {
        display(glPanel, StsCultureObject2D.D2);
    }
    
    public void display(StsGLPanel glPanel)
    {  
        display(glPanel, StsCultureObject2D.D3);
    }
    
    public void display(StsGLPanel glPanel, byte dimensions)
    {
        if(!getIsVisible()) return;
        GL gl = glPanel.getGL();

        if(displayListChanged)
        {
            clearDisplayList(gl);
            displayListChanged = false;
        }
        if(!checkLoadPlanes(gl))
        	return;
 
        boolean useDisplayLists = currentModel.useDisplayLists;
        if (gridDisplayListNum == 0 && useDisplayLists) // build display list
        {
            gridDisplayListNum = gl.glGenLists(1);
            if (gridDisplayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsMultiAttributeVector.display: Failed to allocate a display list");
                return;
            }

            currentModel.getGlPanel3d().setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glNewList(gridDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
            draw2d3d(cultureDisplayableObject, gl, dimensions);
            gl.glEndList();
            currentModel.getGlPanel3d().resetViewShift(gl);
        }
        else if (useDisplayLists) // use existing display list
        {
        	currentModel.getGlPanel3d().setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glCallList(gridDisplayListNum);
            currentModel.getGlPanel3d().resetViewShift(gl);
        }
        else // immediate mode draw
        {
        	clearDisplayList(gl);
            currentModel.getGlPanel3d().setViewShift(gl, StsGraphicParameters.gridShift);
            draw2d3d(cultureDisplayableObject, gl, dimensions);
            currentModel.getGlPanel3d().resetViewShift(gl);
        }        
    }

    public void draw2d3d(StsCultureDisplayable cultureDisplayable, GL gl, byte dimensions)
    {
        if(volumes == null) return;
        if(volumes.length < 1) return;

        if(colorscale == null)
            initializeColorscale();

        if(planes == null)
        	return;
        
        float[][] line = new float[3][3];  // Start, middle and end of arrow
        double interval = 255;
        int bin = 0;       
        
        int nPoints = 0;
        StsColor color;
        float z = currentZ;
        for(int row = 0; row < numRows; row++)
        {
            for(int col = 0; col < numCols; col++)
            {
                if(nPoints%decimation != 0)
                {
                    nPoints++;
                    continue;
                }
                // Determine center point
                line[1][0] = volumes[0].getXCoor(col);
                line[1][1] = volumes[0].getYCoor(row);
                if(!isPlanar && (dimensions == StsCultureObject2D.D3))  // Project
                {
                    z = cultureDisplayable.getCultureZ(line[1][0], line[1][1]) +
                        ((StsSurface)cultureDisplayable).getOffset();
                }
                line[1][2] = z;

                // Compute the length based on attribute 1
                float lengthValue = 1;
                if(volumes[LENGTH] != null)
                {
                	// Skip null values
                	if(planes[LENGTH][nPoints] == -1)
                	{
                		nPoints++;
                		continue;
                	}
                    lengthValue = ((StsSeismicVolume)volumes[LENGTH]).computeFloatFromByte(planes[LENGTH][nPoints]);
                    if(applyThreshold)
                    {
                        if(lengthValue < lengthThreshold)
                        {
                            nPoints++;
                            continue;
                        }
                    }
                    interval = (volumes[LENGTH].dataMax- volumes[LENGTH].dataMin)/6.0f;
                    lengthValue = (int) ((lengthValue - volumes[LENGTH].dataMin) / interval);
                    lengthValue = lengthValue * scale;
                }

                // Compute the direction based on attribute 2
                float azimuthValue = 0;
                if(volumes[AZIMUTH] != null)
                {
                    azimuthValue = ((StsSeismicVolume)volumes[AZIMUTH]).computeFloatFromByte(planes[AZIMUTH][nPoints]);
                    if(normalizeAzimuth)
                    {
                        interval = (volumes[AZIMUTH].dataMax- volumes[AZIMUTH].dataMin)/360.0f;
                        if((interval > 1.0) || (interval < 0.975))
                        {
                            bin = (int) ((azimuthValue - volumes[AZIMUTH].dataMin) / interval);
                            azimuthValue = bin % 360.0f;
                        }
                    }
                }
                // Compute the color based on attribute 3
                color = new StsColor(StsColor.RED);
                if(volumes[COLOR] != null)
                {                	
                    float colorValue = ((StsSeismicVolume)volumes[COLOR]).computeFloatFromByte(planes[COLOR][nPoints]);
                    interval = (volumes[COLOR].dataMax- volumes[COLOR].dataMin)/(double)(colorscale.getNColors()-1);
                    bin = (int)((colorValue - volumes[COLOR].dataMin)/interval);
                    int colorIdx = bin % colorscale.getNColors();
                    color = colorscale.getStsColor(colorIdx);
                }

                // Compute the end points using the length and azimuth values.
                double dx, mdx;
                double dy, mdy;
                dx = Math.sin(azimuthValue) * lengthValue/2.0f;
                dy = Math.cos(azimuthValue) * lengthValue/2.0f;
                mdx = Math.sin(90+azimuthValue) * lengthValue/6.0f;
                mdy = Math.cos(90+azimuthValue) * lengthValue/6.0f;

                line[0][0] = line[1][0] - (float)dx;
                line[0][1] = line[1][1] - (float)dy;
                line[0][2] = line[1][2];
                line[2][0] = line[1][0] + (float)dx;
                line[2][1] = line[1][1] + (float)dy;
                line[2][2] = line[1][2];
                //System.out.println("Azimuth=" + azimuthValue + " Length=" + lengthValue);

                // Draw the line
                float[][] rectangle = new float[5][3];
                float[][] diamond = new float[5][3];
                float[] ctrPoint = new float[3];
                boolean endPoints = false;

                ctrPoint = line[1];

                StsPoint[] points = getPoints(line);

                switch(symbolType)
                {
                    case LINE:
                        break;
                    case LINEENDS:
                        endPoints = true;
                        break;
                    case DIAMOND:
                        System.arraycopy(line[0],0,diamond[0],0,3);
                        diamond[1][0] = line[1][0] - (float)mdx;
                        diamond[1][1] = line[1][1] - (float)mdy;
                        diamond[1][2] = line[0][2];
                        System.arraycopy(line[2],0,diamond[2],0,3);
                        diamond[3][0] = line[1][0] + (float)mdx;
                        diamond[3][1] = line[1][1] + (float)mdy;
                        diamond[3][2] = line[0][2];
                        System.arraycopy(line[0],0,diamond[4],0,3);
                        points = getPoints(diamond);
                        line = diamond;
                        break;
                    case BOWTIE:
                        rectangle[0][0] = line[2][0] + (float)mdx;
                        rectangle[0][1] = line[2][1] + (float)mdy;
                        rectangle[0][2] = line[0][2];
                        rectangle[1][0] = line[2][0] - (float)mdx;
                        rectangle[1][1] = line[2][1] - (float)mdy;
                        rectangle[1][2] = line[0][2];
                        rectangle[2][0] = line[0][0] + (float)mdx;
                        rectangle[2][1] = line[0][1] + (float)mdy;
                        rectangle[2][2] = line[0][2];
                        rectangle[3][0] = line[0][0] - (float)mdx;
                        rectangle[3][1] = line[0][1] - (float)mdy;
                        rectangle[3][2] = line[0][2];
                        rectangle[4][0] = line[2][0] + (float)mdx;
                        rectangle[4][1] = line[2][1] + (float)mdy;
                        rectangle[4][2] = line[0][2];
                        points = getPoints(rectangle);
                        line = rectangle;
                        break;
                }
                if(dimensions == StsCultureObject2D.D3)
                {
                    if (stroke == SOLID)
                        StsGLDraw.drawLineStrip(gl, color, line, 2);
                    else
                        StsGLDraw.drawDottedLineStrip(gl, color, line, 2);
                    if(endPoints)
                    {
                        StsGLDraw.drawPoint(gl, line[0], StsColor.DARK_GRAY, 3);
                        StsGLDraw.drawPoint(gl, line[2], StsColor.DARK_GRAY, 3);
                    }
                    if(drawMidpoints) StsGLDraw.drawPoint(gl, ctrPoint, StsColor.DARK_GRAY, 3);
                }
                else
                {
                    if(stroke == SOLID)
                        StsGLDraw.drawLineStrip2d(gl, color, points,  2, 2);
                    else
                        StsGLDraw.drawDottedLineStrip2d(gl, color, points, 2, 2);
                    if(endPoints)
                    {
                        StsGLDraw.drawPoint2d(line[0], StsColor.DARK_GRAY, gl, 3);
                        StsGLDraw.drawPoint2d(line[2], StsColor.DARK_GRAY, gl, 3);
                    }
                    if(drawMidpoints) StsGLDraw.drawPoint2d(new StsPoint(ctrPoint), StsColor.DARK_GRAY, gl, 3);
                }
                nPoints++;
            }
        }
    }

    public StsPoint[] getPoints(float[][] floats)
    {
        StsPoint[] points = new StsPoint[floats.length];
        for(int i=0; i<floats.length; i++)
            points[i] = new StsPoint(floats[i][0],floats[i][1],floats[i][2]);
        return points;
    }
    public int getStroke() { return stroke; }
    public void setStroke(int stroke) { this.stroke = stroke; }
    public boolean getDrawMidpoints() { return drawMidpoints; }
    public void setDrawMidpoints(boolean midpts) 
    { 
    	this.drawMidpoints = midpts; 
    	clearDisplayList();
        currentModel.viewObjectRepaint(this, this);
    }
    public String[] getTypeStrings() { return typeStrings; }

    private void initializePanel()
    {
        StsObject[] modelSurfaces = (StsObject[])currentModel.getObjectList(StsModelSurface.class);
        StsObject[] surfaces = (StsObject[])currentModel.getObjectList(StsSurface.class);
        int nModelSurfaces = modelSurfaces.length;
        int nSurfaces = surfaces.length;
        Object[] listItems = new Object[nModelSurfaces + nSurfaces + 1];
        listItems[0] = currentModel.getProject();
        System.arraycopy(surfaces, 0, listItems, 1, nSurfaces);
        System.arraycopy(modelSurfaces, 0, listItems, nSurfaces+1, nModelSurfaces);
        cultureDisplayableList.setListItems(listItems);
    }

    public StsFieldBean[] getDisplayFields()
    {
        StsComboBoxFieldBean symbolListBean = new StsComboBoxFieldBean(StsMultiAttributeVector.class, "symbolString", "Symbol:", typeStrings);
        StsFloatFieldBean thresholdBean = new StsFloatFieldBean(StsMultiAttributeVector.class, "lengthThreshold", volumes[LENGTH].getDataMin(),
            volumes[LENGTH].getDataMax(),"Threshold:", true);
        thresholdBean.setRangeFixStep(volumes[LENGTH].getDataMin(), volumes[LENGTH].getDataMax(),
                                      (volumes[LENGTH].getDataMax()-volumes[LENGTH].getDataMin())/1000.0f);
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsMultiAttributeVector.class, "name", "Name"),
            new StsBooleanFieldBean(StsMultiAttributeVector.class, "isVisible", "Enable"),
            new StsIntFieldBean(StsMultiAttributeVector.class, "scale", 1, 100, "Length Scale:", true),
            new StsIntFieldBean(StsMultiAttributeVector.class, "decimation", 1, 100, "Decimation:", true),
            new StsBooleanFieldBean(StsMultiAttributeVector.class, "applyThreshold", "Apply Threshold"),
            new StsBooleanFieldBean(StsMultiAttributeVector.class, "drawMidpoints", "Draw Midpoints"),
            thresholdBean,
            symbolListBean,
            colorscaleBean
//            new StsIntFieldBean(StsCultureObjectSet2D.class, "width", true, "Line Width:", true)
        };
        return displayFields;
    }

    public String getSymbolString()
    {
        return typeStrings[symbolType];
    }

    public void setSymbolString(String symbolString)
    {
        for(int i=0; i<typeStrings.length; i++)
        {
            if (symbolString.equals(typeStrings[i]))
            {
                symbolType = (byte)i;
                break;
            }
        }
        clearDisplayList();
        dbFieldChanged("symbolType", symbolType);
        currentModel.viewObjectRepaint(this, this);
        return;
    }
    public float getLengthThreshold() { return lengthThreshold; }
    public void setLengthThreshold(float threshold)
    {
        this.lengthThreshold = threshold;
        clearDisplayList();
        dbFieldChanged("lengthThreshold", lengthThreshold);
        currentModel.viewObjectRepaint(this, this);
    }
    public boolean getApplyThreshold() { return applyThreshold; }
    public void setApplyThreshold(boolean applyIt)
    {
        this.applyThreshold = applyIt;
        clearDisplayList();
        dbFieldChanged("applyThreshold", applyThreshold);
        currentModel.viewObjectRepaint(this, this);
    }
    public int getDecimation() { return decimation; }
    public void setDecimation(int decimation)
    {
        this.decimation = decimation;
        clearDisplayList();
        dbFieldChanged("decimation", decimation);
        currentModel.viewObjectRepaint(this, this);
    }
    public int getScale() { return scale; }
    public void setScale(int scale)
    {
        this.scale = scale;
        clearDisplayList();
        dbFieldChanged("scale", scale);
        currentModel.viewObjectRepaint(this, this);
    }
    public int getWidth() { return width; }
    public void setWidth(int width)
    {
        this.width = width;
        dbFieldChanged("width", width);
        currentModel.viewObjectRepaint(this, this);
    }
    public boolean anyDependencies()
    {
        return false;
    }
	
    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields != null) return propertyFields;
        cultureDisplayableList = new StsComboBoxFieldBean(StsMultiAttributeVector.class, "cultureDisplayable", "Map Vectors To:");
        propertyFields = new StsFieldBean[]
        {
            cultureDisplayableList
        };
        return propertyFields;
    }
    public void setCultureDisplayable(StsCultureDisplayable cultureDisplayable)
    {
        cultureDisplayableObject = cultureDisplayable;
        dbFieldChanged("cultureDisplayableObject", cultureDisplayableObject);
        currentModel.viewObjectRepaint(this, this);
    }
    public StsCultureDisplayable getCultureDisplayable() { return cultureDisplayableObject; }
    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        initializePanel();
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        getMultiAttributeVectorClass().selected(this);
    }
    public Object[] getChildren()
    {
        return new Object[0];
    }

}
