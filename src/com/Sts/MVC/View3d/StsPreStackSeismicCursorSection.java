package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.nio.*;
//TODO Need to finally remove StsPreStackCursorSection and create 3 different PreStackVolume subclasses all of which would
//TODO be displayed on StsSeismicCursorSection instead.
public class StsPreStackSeismicCursorSection extends StsCursor3dTexture
{
    /** 3d lineSet currently being displayed */
    public StsPreStackLineSet3d lineSet3d;
    /** data currently being drawn on this cursorSection */
    transient public ByteBuffer planeData = null;
    /** display type on this cursorSection: velocity, stack, semblance, attribute, or none (see StsPreStackLineSetClass.displayModeStrings. */
    transient public byte displayType = DISPLAY_NONE;
    /** plane currently being drawn (recomputed on every displayTexture call */
    transient int nPlane = -1;
    /** Display lists should be used (controlled by View:Display Options) */
    // transient boolean useDisplayLists;
    /** Display lists currently being used for surface geometry */
    // transient boolean usingDisplayLists = false;
    /** indicates whether cursorSection is being drawn in time or depth. */
    transient byte zDomain = StsParameters.TD_NONE;
    /** indicates cursor is currently being displayed in 3d */
    transient boolean is3d = true;
    /** type of stack on this section: none, neighbor, line, or volume */
    transient byte stackOption;

    transient public boolean cropChanged = false;
    transient boolean isPixelMode;
    transient boolean isVisible = true;

    transient StsPreStackLineSetClass lineSetClass;

    static final byte STACK_NONE = StsPreStackLineSet3dClass.STACK_NONE;
    static final byte STACK_NEIGHBORS = StsPreStackLineSet3dClass.STACK_NEIGHBORS;
    static final byte STACK_LINES = StsPreStackLineSet3dClass.STACK_LINES;
    static final byte STACK_VOLUME = StsPreStackLineSet3dClass.STACK_VOLUME;

    static final byte DISPLAY_VELOCITY = StsPreStackLineSetClass.DISPLAY_VELOCITY;
    static final byte DISPLAY_STACKED = StsPreStackLineSetClass.DISPLAY_STACKED;
    static final byte DISPLAY_SEMBLANCE = StsPreStackLineSetClass.DISPLAY_SEMBLANCE;
    static final byte DISPLAY_ATTRIBUTE = StsPreStackLineSetClass.DISPLAY_ATTRIBUTE;
    static final byte DISPLAY_NONE = StsPreStackLineSetClass.DISPLAY_NONE;

    transient private StsTimer timer = new StsTimer();
    transient private boolean runTimer = false;

    /** minimum size in pixels of a crossplot point */
    static private final int minPointSize = 4;

    public StsPreStackSeismicCursorSection()
    {
    }

    public StsPreStackSeismicCursorSection(StsPreStackLineSet3d seismicVolume, StsModel model, StsCursor3d cursor3d, int dir)
    {
        this.lineSet3d = seismicVolume;
        initialize(model, cursor3d, dir);
    }

    public boolean initialize(StsModel model, StsCursor3d cursor3d, int dir)
    {
        super.initialize(model, cursor3d, dir);
        if (lineSet3d == null) return false;
        lineSetClass = (StsPreStackLineSet3dClass) model.getStsClass(StsPreStackLineSet3d.class);
        isPixelMode = getIsPixelMode();
        this.nTextureRows = lineSet3d.maxNTracesPerGather;
        this.nTextureCols = lineSet3d.nSlices;
        lineSet3d.setDirCoordinate(dirNo, dirCoordinate, cursor3d.window);
        return true;
    }

    protected boolean getIsPixelMode()
    {
        return lineSetClass.getIsPixelMode();
    }

    /*
        private void initShader()
        {
                boolean useShader = seismicClass.getContourColors() && StsJOGLShader.canUseShader;
                if(useShader)
                    shader = StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
                else
                    shader = StsJOGLShader.NONE;
        }

        private boolean shaderChanged()
        {
            boolean usingShader = textureTiles.shader != StsJOGLShader.NONE;
            boolean useShader = seismicClass.getContourColors() && StsJOGLShader.canUseShader;
            if(useShader == usingShader) return false;
            textureChanged = true;
            if(useShader)
                shader = StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
            else
                shader = StsJOGLShader.NONE;
            return true;
        }
    */
    public int getDefaultShader()
    { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }

    public boolean getUseShader() { return lineSetClass.getContourColors(); }

    /**
     * This puts texture display on delete list.  Operation is performed
     * at beginning of next draw operation.
     */
    public boolean textureChanged()
    {
        // TODO:  if this is a zSlice and we are displaying attribute, the texture hasn't changed
        // TODO: so we should check this and return with setting textureChanged
        if (lineSet3d == null) return false;
        if (textureTiles != null)
        {
            if (debug) StsException.systemDebug(this, "textureChanged", debugMessage("deleting textureTileSurface"));
            textureChanged = true;
            planeData = null;
        }
        return textureChanged;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void cropChanged()
    {
//		if(textureTiles == null && !checkTextureTiles())return;
        if (textureTiles == null) return;
        textureChanged();
        textureTiles.cropChanged();
    }

    public void subVolumeChanged()
    {
    }

    public Class getDisplayableClass() { return StsPreStackLineSet3d.class; }

    public boolean canDisplayClass(Class c) { return StsPreStackLineSet3d.class.isAssignableFrom(c); }

    public boolean isDisplayableObject(Object object)
    {
        return (object instanceof StsPreStackLineSet3d);
    }

    public boolean isDisplayingObject(Object object)
    {
        if (lineSet3d == object) return true;
        if (this.textureTiles != null && object == textureTiles.cropVolume) return true;
        return false;
    }

    /*
      public boolean setDirCoordinate(int dirNo, float dirCoordinate)
      {
       if(line2d.setDirCoordinate(dirNo, dirCoordinate)) return true;
       clearTextureTileSurface(glPanel3d);
 //		if (textureTiles == null) checkTextureTiles();
 //		textureTiles.setTilesDirCoordinate(dirCoordinate);
       return true;
      }
      */

    public boolean setDirCoordinate(float dirCoordinate)
    {
        if (lineSet3d != null) lineSet3d.setDirCoordinate(dirNo, dirCoordinate, cursor3d.window);
        if (this.dirCoordinate != dirCoordinate)
        {
            textureChanged();
            unlockTexturePlanes();
            if(debug && this.dirCoordinate != dirCoordinate) StsException.systemDebug(this, "setDirCoordinate", "dir " + dirNo + " coor " + this.dirCoordinate + " changed to " + dirCoordinate);
            this.dirCoordinate = dirCoordinate;
            geometryChanged();
        }
        return true;
    }

    private void unlockTexturePlanes()
    {
        lineSet3d.unlockPlanes(dirNo, dirCoordinate);
    }

    public boolean setObject(Object object)
    {
        if(object == lineSet3d) return false;
        if (isDisplayableObject(object))
        {
            lineSet3d = (StsPreStackLineSet3d) object;
            lineSet3d.setDirCoordinate(dirNo, dirCoordinate, cursor3d.window);
            isVisible = true;
        }
        else if (object != null)
            return false;
        else // object == null
            isVisible = false;
        textureChanged();
        return true;
    }

    public Object getObject()
    {
        return lineSet3d;
    }

    public void displayTexture(StsGLPanel3d glPanel3d, boolean is3d, StsCursorSection cursorSection)
    {
//		if(!line2d.canDisplayZDomain())return;
        if (lineSet3d == null) return;
        if (!isVisible) return;
        if (glPanel3d.getCursor3d().getIsDragging(dirNo) && (dirNo != StsCursor3d.ZDIR)) return;
        if (!lineSet3d.getIsVisibleOnCursor()) return;

        GL gl = glPanel3d.getGL();

        if(debug)
        {
            int nPlaneOld = nPlane;
            nPlane = lineSet3d.getCursorPlaneIndex(dirNo, dirCoordinate);
            if(nPlaneOld != nPlane) StsException.systemDebug(this, "displayTexture", "dir " + dirNo + " nPlane " + nPlaneOld + " changed to " + nPlane);
        }
        else
            nPlane = lineSet3d.getCursorPlaneIndex(dirNo, dirCoordinate);

        checkTextureAndGeometryChanges();
        if (!initializeTextureTiles(glPanel3d)) return;

        if (this.is3d != is3d)
        {
            //if (debug)
//                StsException.systemDebug(this, "displayTexture", "Deleting seismic texture. " + textureTiles + " texture changed " + textureChanged);
//            textureChanged = true;
            geometryChanged = true;
            this.is3d = is3d;

            if (!is3d)
            {
                StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getView();
                textureTiles.setAxesFlipped(viewCursor.axesFlipped);
            }
        }
         // if geometryChanged, delete the displayLists and reset flag. Clear textures has they must be refetched or recomputed
        // Missing displayLists are a flag themselves indicated they need to be rebuilt if being used.
        if(geometryChanged)
        {
            textureChanged = true;
            deleteDisplayLists(gl);
            geometryChanged = false;
        }
        // if texturesChanged, clear the textureTiles but don't reset flag as it will be used to initiate recomputing textures
        if (textureChanged && displayType != DISPLAY_ATTRIBUTE)
             textureTiles.deleteTextures(gl);



        /** plane for this direction and coordinate can't be drawn because it is cropped out */
        if (textureTiles.isDirCoordinateCropped()) return;

        textureTiles.checkBuildDisplayLists(gl, is3d);

        if (textureTiles.isDirCoordinateCropped()) return;

        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_BLEND);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);

            if (!checkSetPlaneData(gl))
            {
                textureChanged = false;
                return;
            }
            if (cursorSection.subVolumePlane != null)
            {
                // planeData = applySubVolume(planeData, subVolumePlane);
            }

            if (debug) StsException.systemDebug(this, "displayTexture", debugMessage("textureChanged: " + textureChanged));
            if (is3d)
            {
                if (textureChanged)
                {
                    textureTiles.displayTiles(this, gl, isPixelMode, planeData, nullByte);
                    textureChanged = false;
                }
                else
                    textureTiles.displayTiles(this, gl, isPixelMode, (ByteBuffer) null, nullByte);
            }
            else
            {
                StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getView();
                if (textureChanged)
                {
//				System.out.println("Displaying changed texture 2d.");
                    textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, planeData, nullByte);
                    textureChanged = false;
                }
                else
                {
//				System.out.println("Displaying old texture 2d.");
                    textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, (ByteBuffer) null, nullByte);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsPreStackSeismicCursorSection.displayTexture() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glDisable(GL.GL_BLEND);
            gl.glEnable(GL.GL_LIGHTING);
            if (textureTiles.shader != StsJOGLShader.NONE)
                StsJOGLShader.disableARBShader(gl);
        }
    }

    private void checkTextureAndGeometryChanges()
    {
        if (isPixelMode != lineSetClass.getIsPixelMode())
        {
            textureChanged = true;
            isPixelMode = !isPixelMode;
        }
        if (lineSetChanged())
        {
            textureChanged = true;
        }
        if (stackOptionChanged())
            textureChanged = true;

        if (displayTypeChanged())
            textureChanged = true;

        byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
        if (projectZDomain != zDomain)
        {
            geometryChanged = true;
            zDomain = projectZDomain;
        }
    }

    protected boolean initializeTextureTiles(StsGLPanel3d glPanel3d)
    {
        // if (!glPanel3d.initialized) return false;
        if (lineSet3d == null) return false;

        if (textureTiles == null)
        {
            StsCropVolume subVolume = model.getProject().getCropVolume();
            textureTiles = StsTextureTiles.constructor(model, this, dirNo, lineSet3d, isPixelMode, subVolume);
            if (textureTiles == null) return false;
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
            geometryChanged = true;
        }
        else if (!textureTiles.isSameSize(lineSet3d))
        {
            textureTiles.constructTiles(lineSet3d);
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
            geometryChanged = true;
        }
        else if (textureTiles.shaderChanged())
        {
            textureChanged = true;
        }
        textureTiles.setTilesDirCoordinate(dirCoordinate);
        return true;
    }

    /**
     * We can display these types of data on all 3 planes: VELOCITY, STACK, and SEMBLANCE.
     * On the z-plane we can also display ATTRIBUTE data.
     * So on the X and Y planes, we attempt to display the type selected (VELOCITY, STACK, or SEMBLANCE).
     * On the Z plane, the data displayed is selected with this priority:
     * 1. If STACK, SEMBLANCE, or VELOCITY is available (only possible if full volume has been computed), display it.
     * 2. If not available and ATTRIBUTE is available, display it.
     * 3. Display NONE.
     *
     * @param gl
     * @return
     */
    private boolean checkSetPlaneData(GL gl)
    {
        if(displayType == DISPLAY_NONE) return false;
        boolean dataOk = false;
        if (displayType == DISPLAY_STACKED)
            dataOk = checkSetStackPlaneData(nPlane);
        else if (displayType == DISPLAY_SEMBLANCE)
            dataOk = checkSetSemblancePlaneData(nPlane);
        else if (displayType == DISPLAY_VELOCITY)
            dataOk = checkSetVelocityPlaneData(nPlane);
        else if(displayType == DISPLAY_ATTRIBUTE)
            dataOk = checkSetAttributePlaneData(nPlane);
        if(!dataOk) return false;
        return lineSet3d.setGLColorList(gl, false, displayType, textureTiles.shader);
    }

    private boolean checkSetVelocityPlaneData(int nPlane)
    {
        if (dirNo == ZDIR) return false;
        if (textureChanged)
        {
            textureTiles.setTilesDirCoordinate(dirCoordinate);
            planeData = lineSet3d.computeVelocityPlane(dirNo, nPlane);
        }
        return planeData != null;
    }

    private boolean checkSetStackPlaneData(int nPlane)
    {
        StsPreStackVolume stackedVolume = lineSet3d.getStackedVolume();
        if (!textureChanged)
        {
            if (lineSet3d.velocityModel.interpolation.upToDate)
                textureChanged = !stackedVolume.isPlaneOK(dirNo, nPlane);
        }
        if (textureChanged)
        {
            textureTiles.setTilesDirCoordinate(dirCoordinate);
            planeData = lineSet3d.computeStackedPlane(dirNo, nPlane);
        }
        return planeData != null && stackedVolume.isPlaneOK(dirNo, nPlane);
    }

    private boolean checkSetSemblancePlaneData(int nPlane)
    {
        StsPreStackVolume semblanceVolume = lineSet3d.getStackedVolume();
        if (!textureChanged)
        {
            if (lineSet3d.velocityModel.interpolation.upToDate)
                textureChanged = !semblanceVolume.isPlaneOK(dirNo, nPlane);
        }
        if (textureChanged)
        {
            textureTiles.setTilesDirCoordinate(dirCoordinate);
            planeData = lineSet3d.computeSemblancePlane(dirNo, nPlane);
        }
        return planeData != null && semblanceVolume.isPlaneOK(dirNo, nPlane);
    }

    private boolean checkSetAttributePlaneData(int nPlane)
    {
        textureTiles.setTilesDirCoordinate(dirCoordinate);
        textureChanged = lineSet3d.isAttributeTextureChanged();
        planeData = lineSet3d.getAttributeByteBuffer(dirNo);
        return planeData != null;
    }

    /*
        private boolean displayAttributeChanged()
        {
            String currentDisplayAttribute = lineSet.getCurrentDisplayAttribute(dirNo);
            if(currentDisplayAttribute == displayAttribute) return false;
            displayAttribute = currentDisplayAttribute;
            return true;
        }
    */
    private boolean stackOptionChanged()
    {
        byte currentStackOption = lineSetClass.getStackOption();
        if (currentStackOption == stackOption) return false;
        stackOption = currentStackOption;
        return true;
    }

   /** Prestack volumes (velocity, stack, semblance) can always be displayed on X or Y (vertical) planes.
    *  If they haven't been computed for these planes, they will be computed in this display pass.
    *  In this case, we need only check that the dirCoordinate is legitimate.
    *  Stack and semblance can't be displayed on Z (horizontal) plane unless volumes are already computed.
    *  This can be checked by the status of the planes: volume.isPlaneOk(dir, coor).
    *  If not ok, the displayType will be set to ATTRIBUTE if it is a z plane, otherwise set to NONE.
    */
    private boolean displayTypeChanged()
    {
        byte currentDisplayType = lineSetClass.getDisplayType();
        boolean dataOk = lineSet3d.isVolumePlaneOk(currentDisplayType, dirNo, dirCoordinate);
        if (!dataOk)
        {
            boolean hasAttributeData = dirNo == StsCursor3d.ZDIR && lineSet3d.hasAttributeData();
            if (hasAttributeData)
                currentDisplayType = DISPLAY_ATTRIBUTE;
            else
                currentDisplayType = DISPLAY_NONE;
        }

        if (displayType != currentDisplayType)
        {
            displayType = currentDisplayType;
            return true;
        }
        else
            return false;
    }

    private boolean lineSetChanged()
    {
        StsPreStackLineSet currentLineSet = lineSetClass.getCurrentProjectLineSet();
        if (currentLineSet == null) return true;
        if (!(currentLineSet instanceof StsPreStackLineSet3d))
        {
            lineSet3d = null;
            return true;
        }
        if (lineSet3d == currentLineSet) return false;
        lineSet3d = (StsPreStackLineSet3d) currentLineSet;
        return true;
    }

    /**
     * The subVolumePlane has the dimensions of the cursorBoundingBox it cuts through;
     * The planeData has the dimensions of the corresponding seismic/virtual volume.
     * If they are the same size, apply the planeData directly to the subVolumePlane;
     * if not get the grid coordinates of the planeData relative to the subVolume plane
     * and apply.
     */
    /*
      protected byte[] applySubVolume(byte[] planeData, byte[] subVolumePlane)
      {
       int nCursorCols, nCursorRows, subVolRowStart, subVolRowEnd, subVolColStart, subVolColEnd;
       int subVolIndex, cursorIndex;
       int cursorRow, cursorCol;
       int subVolRow, subVolCol;

       if (planeData == null)
       {
        return null;
       }
       if (subVolumePlane == null)
       {
        return planeData;
       }
       if (!line2d.getPreStackSeismicVolumeClass().getDisplayOnSubVolumes())
       {
        return planeData;
       }

       byte[] subVolumeData = new byte[subVolumePlane.length];

       try
       {

 //        if (line2d.congruentWith(cursorBoundingBox))
        {
         if (planeData.length == subVolumePlane.length)
         {
          for (int n = 0; n < planeData.length; n++)
          {
           if (subVolumePlane[n] == 0)
           {
            subVolumeData[n] = -1;
           }
           else
           {
            subVolumeData[n] = planeData[n];
           }
          }
         }
         else
         {
          StsSubVolumeClass subVolumeClass = (StsSubVolumeClass) glPanel3d.model.getStsClass(StsSubVolume.class);
          StsGridBoundingBox subVolumeBoundingBox = subVolumeClass.getBoundingBox();
          nCursorCols = cursorBoundingBox.getNCursorCols(dirNo);
          nCursorRows = cursorBoundingBox.getNCursorRows(dirNo);
          subVolRowStart = subVolumeBoundingBox.getCursorRowMin(dirNo);
          subVolRowEnd = subVolumeBoundingBox.getCursorRowMax(dirNo);
          subVolColStart = subVolumeBoundingBox.getCursorColMin(dirNo);
          subVolColEnd = subVolumeBoundingBox.getCursorColMax(dirNo);
          subVolIndex = 0;
          cursorIndex = 0;
          // make rows before subVolRowStart transparent
          for (cursorRow = 0; cursorRow < subVolRowStart; cursorRow++)
          {
           for (cursorCol = 0; cursorCol < nCursorCols; cursorCol++, cursorIndex++)
           {
            subVolumeData[cursorIndex] = -1;
           }
          }
          for (subVolRow = subVolRowStart; subVolRow <= subVolRowEnd; subVolRow++)
          {
           cursorIndex = subVolRow * nCursorCols + subVolColStart;
           for (subVolCol = subVolColStart; subVolCol <= subVolColEnd; subVolCol++, subVolIndex++,
             cursorIndex++)
           {
            if (subVolumePlane[subVolIndex] == 0)
            {
             subVolumeData[cursorIndex] = -1;
            }
            else
            {
             subVolumeData[cursorIndex] = planeData[cursorIndex];
            }
           }
          }
          // make rows after subVolRowStart transparent
          for (cursorRow = subVolRowEnd + 1; cursorRow < nCursorRows; cursorRow++)
          {
           for (cursorCol = 0; cursorCol < nCursorCols; cursorCol++, cursorIndex++)
           {
            subVolumeData[cursorIndex] = -1;
           }
          }
         }
         return subVolumeData;
        }
       }
       catch (Exception e)
       {
        StsException.outputException("StsPreStackSeismicCursorSection.applySubVolume() failed.", e, StsException.WARNING);
        return null;
       }
      }
      */
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        lineSet3d.drawTextureTileSurface(tile, gl, dirNo, is3d);
    }

    /*
      public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
      {
       line2d.drawTextureTileSurface(tile, gl, dirNo, false);
      }
      */
    public void display(StsGLPanel3d glPanel3d, boolean is3d)
    {
        if (lineSet3d == null) return;
        if (!lineSet3d.canDisplayZDomain()) return;

        GL gl = glPanel3d.getGL();
        if (is3d)
        {
//			displayCrossplotPoints(gl);
        }
        else
        {
            StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getView();
            GLU glu = glPanel3d.getGLU();
            displayWiggleTraces(glPanel3d, gl, viewCursor.axisRanges);
        }
    }

    private void displayWiggleTraces(StsGLPanel3d glPanel3d, GL gl, float[][] axisRanges)
    {
        if (!displayWiggles()) return;
        if (planeData == null) return;

        int nSlices = lineSet3d.nSlices;
        int nTraces = 0;
        float horizMin = 1.0f;
        float height = cursorBoundingBox.zMax - cursorBoundingBox.zMin;
        float verticalMin = cursorBoundingBox.zMin;
        float horizInc = 1.0f;

        // If density is less than 1:4 traces to pixels, display wiggles
        if (((axisRanges[0][1] - axisRanges[0][0]) / horizInc) > glPanel3d.getWidth() / getWiggleToPixelRatio())
            return;

        if (planeData == null) return;

        gl.glPushMatrix();
        gl.glLoadIdentity();
        float horizScale = 3 * horizInc / 254;
        float verticalScale = height / (nSlices - 1);
        horizInc /= horizScale;
        gl.glScalef(horizScale, verticalScale, 1.0f);
        gl.glTranslatef(horizMin / horizScale, verticalMin / verticalScale, 0.0f);
        StsColor.BLACK.setGLColor(gl);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glLineWidth(0.5f);
        boolean isDrawing = false;
        int n = 0;

        planeData.rewind();
        for (int t = 0; t < nTraces; t++)
        {
            for (int s = 0; s < nSlices; s++, n++)
            {
                byte value = planeData.get();
                if (value == -1)
                {
                    if (isDrawing)
                    {
                        gl.glEnd();
                        isDrawing = false;
                    }
                    continue;
                }
                else if (!isDrawing)
                {
                    isDrawing = true;
                    gl.glBegin(GL.GL_LINE_STRIP);
                }
                gl.glVertex2i(StsMath.unsignedByteToSignedInt(value), s);
            }
            if (isDrawing) gl.glEnd();
            isDrawing = false;
            gl.glTranslatef(horizInc, 0.0f, 0.0f);
        }
        gl.glPopMatrix();
        gl.glEnable(GL.GL_LIGHTING);
    }

    private boolean displayWiggles()
    {
        StsPreStackLineSet3dClass aClass = (StsPreStackLineSet3dClass) model.getStsClass(lineSet3d.getClass());
        return aClass.getDisplayWiggles();
    }

    private int getWiggleToPixelRatio()
    {
        StsPreStackLineSet3dClass aClass = (StsPreStackLineSet3dClass) model.getStsClass(lineSet3d.getClass());
        return aClass.getWiggleToPixelRatio();
    }

    public void displayCrossplotPoints(StsGLPanel3d glPanel3d, GL gl, StsCrossplotPoint[] crossplotPoints)
    {
        if (crossplotPoints == null)
        {
            return;
        }
        int nPoints = crossplotPoints.length;

        gl.glDisable(GL.GL_LIGHTING);

        glPanel3d.setViewShift(gl, 1.0);
        for (int n = 0; n < nPoints; n++)
        {
            StsGLDraw.drawPoint(gl, crossplotPoints[n].getVolumeXYZ(), StsColor.BLACK, minPointSize + 2);

        }
        glPanel3d.setViewShift(gl, 2.0);
        for (int n = 0; n < nPoints; n++)
        {
            StsGLDraw.drawPoint(gl, crossplotPoints[n].getVolumeXYZ(), crossplotPoints[n].stsColor, minPointSize);

        }
        glPanel3d.resetViewShift(gl);

        gl.glEnable(GL.GL_LIGHTING);

        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glPopMatrix();
        gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

    public String propertyReadout(StsPoint point)
    {
        StringBuffer stringBuffer = null;
        StsObject[] volumes = model.getObjectList(lineSet3d.getClass());
        for (int i = 0; i < volumes.length; i++)
        {
            StsPreStackLineSet3d volume = (StsPreStackLineSet3d) volumes[i];
            if (!volume.getReadoutEnabled() && volume != lineSet3d)
            {
                continue;
            }
            if (stringBuffer == null)
            {
                stringBuffer = new StringBuffer();
            }
            stringBuffer.append(" " + volume.getName() + ": " +
                volume.getScaledValue((byte) volume.getPlaneValue(point.v)));
        }
        if (stringBuffer == null)
        {
            return null;
        }
        return stringBuffer.toString();
    }

    public String logReadout2d(StsPoint2D point)
    {
        int[] dataRowCol = getDataRowCol(point);
        return logReadout2d(dataRowCol);
    }

    public String logReadout2d(float[] rowColF)
    {
        int row = Math.round(rowColF[0]);
        int col = Math.round(rowColF[1]);
        int[] rowCol = new int[]
            {row, col};
        return logReadout2d(rowCol);
    }

    public String logReadout2d(int[] rowCol)
    {
        String valueString;

        if (rowCol != null && isInRange(rowCol[0], rowCol[1]))
        {
            if (planeData == null)
            {
                return "no data";
            }
            byte byteValue = planeData.get(rowCol[0] * nTextureCols + rowCol[1]);
            float value = lineSet3d.getScaledValue(byteValue);
            return Float.toString(value);
        }
        else
        {
            return "not in range";
        }
        /*
            int row = rowCol[0];
            int col = rowCol[1];
            if(dirNo == StsCursor3d.XDIR)
            {
             float yLabel = line2d.getNumFromIndex(StsCursor3d.YDIR, (float)col);
             float zLabel = line2d.getNumFromIndex(StsCursor3d.ZDIR, (float)row);
             return new String("Line: " + yLabel + " Slice: " + zLabel + " value: " + valueString);
            }
            else if(dirNo == StsCursor3d.YDIR)
            {
             float xLabel = line2d.getNumFromIndex(StsCursor3d.XDIR, (float)col);
             float zLabel = line2d.getNumFromIndex(StsCursor3d.ZDIR, (float)row);
             return new String("Crossline: " + xLabel + " Slice: " + zLabel + " Value: " + valueString);
            }
            else
            {
             float xLabel = line2d.getNumFromIndex(StsCursor3d.XDIR, (float)col);
             float yLabel = line2d.getNumFromIndex(StsCursor3d.YDIR, (float)row);
             return new String("Line: " + yLabel + " Crossline: " + xLabel + " Value: " + valueString);
            }
           */
    }

    public int[] getDataRowCol(StsPoint2D point)
    {
        int row, col;
        /*
            if (dirNo == StsCursor3d.XDIR) // first row is vertically down
            {
             row = (int) line2d.getRowCoor(point.x);
             col = (int) line2d.getSliceCoor(point.y);
            }
            else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
            {
             row = (int) line2d.getColCoor(point.x);
             col = (int) line2d.getSliceCoor(point.y);
            }
            else // dirNo == ZDIR  first row is horizontally across
           */
        {
            row = (int) lineSet3d.getRowCoor(point.y);
            col = (int) lineSet3d.getColCoor(point.x);
        }
        if (isInRange(row, col))
        {
            return new int[]
                {row, col};
        }
        else
        {
            return null;
        }
    }

    /** given 2D coordinates on this 2D cursorSection, return the row col coordinates */
    public float[] getDataRowColF(float x, float y)
    {
        float rowF, colF;
        /*
            if (dirNo == StsCursor3d.XDIR) // first row is vertically down
            {
             rowF = line2d.getRowCoor(x);
             colF = line2d.getSliceCoor(y);
            }
            else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
            {
             rowF = line2d.getColCoor(x);
             colF = line2d.getSliceCoor(y);
            }
            else // dirNo == ZDIR  first row is horizontally across
           */
        {
            rowF = lineSet3d.getRowCoor(y);
            colF = lineSet3d.getColCoor(x);
        }
        if (isInRange(rowF, colF))
        {
            return new float[]
                {rowF, colF};
        }
        else
        {
            return null;
        }
    }

    /** given 3D coordinates on this 3d cursorSection, return the row col coordinates */
    public float[] getDataRowColF(StsPoint point)
    {
        float rowF, colF;
        /*
            if (dirNo == StsCursor3d.XDIR) // first row is vertically down
            {
             rowF = line2d.getRowCoor(point.getY());
             colF = line2d.getSliceCoor(point.getZ());
            }
            else if (dirNo == StsCursor3d.YDIR) // first row is vertically down
            {
             rowF = line2d.getColCoor(point.getX());
             colF = line2d.getSliceCoor(point.getZ());
            }
            else // dirNo == ZDIR  first row is horizontally across
           */
        {
            rowF = lineSet3d.getRowCoor(point.getY());
            colF = lineSet3d.getColCoor(point.getX());
        }
        if (isInRange(rowF, colF))
        {
            return new float[]
                {rowF, colF};
        }
        else
        {
            return null;
        }
    }

    public String getName()
    {
        return "Cursor view[" + dirNo + "] of: " + lineSet3d.getName();
    }
}