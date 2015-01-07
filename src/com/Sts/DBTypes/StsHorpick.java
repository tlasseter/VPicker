package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.Horpick.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.nio.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

public class StsHorpick extends StsMainObject implements StsTreeObjectI
{
    /** array of seeds created by StsHorPickAction */
    protected StsObjectRefList patches = null;
    /** colorscale used for correlCoefs */
    protected StsColorscale corCoefsColorscale = null;
    /** color of seed points and resulting surface */
    protected StsColor stsColor = new StsColor(StsColor.GRAY);
    /** line2d seed points are picked on */
    protected StsSeismicVolume seismicVolume;
    /** surface constructed from pick points by StsAutoTrackAction */
    protected StsSurface surface = null;
    /** indicates picking is being run iteratively from autoCorMax to autoCorMin; otherwise run once at minCorrel */
    protected boolean isIterative = false;
    /** manual picking minimum acceptable cross-correlation */
    public float manualCorMin = defaultMinCorrel;
    /** iterative picking operation max correl */
    public float autoCorMax = defaultAutoCorMax;
    /** iterative picking operation max correl */
    public float autoCorMin = defaultAutoCorMin;
    /** iterative picking operation max correl */
    public float autoCorInc = defaultAutoCorInc;
    /** minCorrel filter optionally set (removes correl values below this threshold) */
    public float minCorrelFilter = 0.0f;
    /** applyMinCorrelFilter indicates we want surface correlation values below min turned off */
    public boolean applyMinCorrelFilter = false;
    /** index of patch at each row, col */
    protected byte[] patchIndexes = null;
    /** corCoefs for this surface */
    protected float[][] corCoefs = null;

    /** texture with patch colors */
    transient PatchColorSurfaceTexture patchColorSurfaceTexture;
    /** texture with corCoefs */
    transient CorrelCoefsSurfaceTexture correlCoefsSurfaceTexture;
    /**
     * the current surfaceTexture associated with this horpick being isVisible on this surface;
     * If the defaultTextureDisplayable, then no texture from this horpick is being isVisible and
     * the surface is displaying its own current textureDisplayable.
     */
    transient StsSurfaceTexture surfaceTexture = null;

    /** list of available surfaceTexture objects for horpick surface */
    transient StsSurfaceTexture[] surfaceTextureList = null;

    transient StsHorpickWizard wizard = null;

    static public final String displayPropertyNone = "None";
    static public final String displayPropertyPatchColor = "Patch Colors";
    static public final String displayPropertyCorrelCoefs = "Correlation Coefficients";

    static public final float defaultAutoCorMax = 0.98f;
    static public final float defaultAutoCorMin = 0.90f;
    static public final float defaultAutoCorInc = 0.02f;
    static public final float defaultMinCorrel = 0.95f;

    static StsObjectPanel objectPanel = null;

    static public transient StsFieldBean[] displayFields = null;

    static StsComboBoxFieldBean textureDisplayableListBean;
    static StsEditableColorscaleFieldBean colorscaleBean;
    static public transient StsFieldBean[] propertyFields = null;

    public StsHorpick()
    {
    }

    private StsHorpick(String name, StsColor color, StsHorpickWizard wizard)
    {
        setName(name);
        this.stsColor = color;
        this.wizard = wizard;
        patches = StsObjectRefList.constructor(5, 5, "patches", this);
//		StsPickPatch.initializeNumber();
//        correlCoefsSurfaceTexture = new CorrelCoefsSurfaceTexture(this);
//        createCorCoefsColorscale();
    }

    static public StsHorpick constructor(String name, StsColor color, StsHorpickWizard wizard)
    {
        try
        {
            return new StsHorpick(name, color, wizard);
        }
        catch (Exception e)
        {
            StsException.outputException("StsHorpick.constructor() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsSurfaceTexture getSurfaceTexture()
    {
        initSurfaceTextureList();
        return surfaceTexture;
    }

    public void initSurfaceTextureList()
    {
        if (surfaceTextureList != null) return;

        int nListItems = 3;
        surfaceTextureList = new StsSurfaceTexture[nListItems];
        surfaceTextureList[0] = surface.surfaceDefaultTexture;
        surfaceTextureList[1] = patchColorSurfaceTexture;
        surfaceTextureList[2] = correlCoefsSurfaceTexture;
        textureDisplayableListBean.setListItems(surfaceTextureList);
//        textureDisplayableListbean.setSelectedIndex(0);
    }

    public void setSurfaceTexture(StsSurfaceTexture surfaceTexture)
    {
        if (surfaceTexture == null) return;
        if (this.surfaceTexture != surfaceTexture)
        {
            this.surfaceTexture = surfaceTexture;
            if (surface != null)
            {
                surface.setSurfaceTexture(surfaceTexture, null);
            }
        }
        colorscaleBean.setHistogram(surfaceTexture.getHistogram());
        StsColorscale colorscale = surfaceTexture.getColorscale();
        colorscaleBean.setValueObject(colorscale);
    }

    // doesn't work properly
    private void checkEnableCorrelCoefBean()
    {
        /*
                   if(objectPanel == null) return;
                   boolean correlCoefsEnabled = newTextureDisplayable == correlCoefsTextureDisplayable;
                   if(correlCoefsEnabled)
              objectPanel.addDisplayField(colorscaleBean);
                   else
              objectPanel.removeDisplayField(colorscaleBean);
           */
    }

    public boolean initialize(StsModel model)
    {
        if (corCoefsColorscale != null)
        {
            corCoefsColorscale.initialize();
            float[] range = getCorCoefRange();
            corCoefsColorscale.setRange(range);
//            corCoefsColorscale.setEditRange(minCorrel, maxCorrel);
        }
        if (surface != null)
        {
            initializeSurfaceTextures();
        }
        return true;
    }

    public void setStsColor(StsColor color)
    {
        if(this.stsColor == color) return;
        stsColor = color;
        StsChangeCmd cmd = new StsChangeCmd(this, stsColor, "stsColor", false);
        currentModel.getCreateTransactionAddCmd("colorChange", cmd);
        currentModel.win3dDisplayAll();
    }

    public StsColor getStsColor()
    {
        if (stsColor == null) return null;
        return stsColor;
    }

    public boolean addPatch(StsPickPatch patch)
    {
        patches.add(patch);
        patchesChanged();
        return true;
    }

    private void patchesChanged()
    {
        surface.setTextureChanged(patchColorSurfaceTexture);
        surface.colorListChanged(patchColorSurfaceTexture);
    }

    public boolean deletePatch(StsPickPatch patch)
    {
        patchesChanged();
        removePickPatchGrid(patch, true);
        patches.delete(patch);
        patch.delete();
        return true;
    }

    public void deleteAllPatches()
    {
        patchesChanged();
        patches.deleteAllElements();
    }

    public boolean deletePicks(StsPickPatch patch)
    {
        removePickPatchGrid(patch, false);
//        saveSurface();
        return true;
    }

    public void deleteSurface()
    {
        if (surface == null) return;
        surface.textureChanged();
        surface.deleteGrid();
    }

    public StsObjectRefList getPatches()
    {
        return patches;
    }

    /*
         public void addPatchToSurface(StsPickPatch patch)
         {
             surface.setHasNulls(true);
             if(surface.getPointsNull() == null)
             {
                 surface.setPointsZ(StsMath.copyFloatArray(patch.getPointsZ()));
                 surface.setPointsNull(StsMath.copyByteArray(patch.getPointsNull()));
             }
             else
             {
                 int nRows = surface.getNRows();
                 int nCols = surface.getNCols();
                 float[][] surfacePointsZ = surface.getPointsZ();
                 byte[][] surfacePointsNull = surface.getPointsNull();
                 float[][] patchPointsZ = patch.getPointsZ();
                 byte[][] patchPointsNull = patch.getPointsNull();
                 for(int row = 0; row < nRows; row++)
                 {
                     for(int col = 0; col < nCols; col++)
                     {
                         byte patchPointNull = patchPointsNull[row][col];
                         if(patchPointNull == StsSurface.NOT_NULL)
                         {
 //                        if(multiValueAction == MULTIVALUE_REPLACE || surfacePointsNull[row][col] != StsSurface.NOT_NULL)
 //                        {
                                 surfacePointsNull[row][col] = patchPointsNull[row][col];
                                 surfacePointsZ[row][col] = patchPointsZ[row][col];
 //                        }
                         }
                     }
                 }
             }
 //        surface.interpolateNullZPoints();
             surface.constructGrid();
         }
      */

    private void removePickPatchGrid(StsPickPatch pickPatch, boolean removeSeed)
    {

        if (patches.getSize() == 1)
        {
            deleteArrays();
            deleteSurface();
        }
        else
        {
            if (patchIndexes == null)
            {
                return;
            }
            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            float[][] surfacePointsZ = surface.getPointsZ();
            byte[][] surfacePointsNull = surface.getPointsNull();
            int pickPatchIndex = pickPatch.getIndex();
            int n = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++, n++)
                {
                    if (patchIndexes[n] == -1)
                    {
                        continue;
                    }
                    int index = StsMath.signedByteToUnsignedInt(patchIndexes[n]);
                    if (index == pickPatchIndex)
                    {
                        surfacePointsNull[row][col] = StsSurface.SURF_GAP_NOT_FILLED;
//                        surfacePointsZ[row][col] = 0.0f;
                        patchIndexes[n] = -1;
                        corCoefs[row][col] = 0.0f;
                        surface.checkClearNeighbors(row, col);
                    }
                    /*
                         else if (removeSeed && index > pickPatchIndex)
                         {
                             patchIndexes[n] = (byte) (index - 1);
                         }
                     */
                }
            }
//            surface.setTextureChanged(true);

//            surface.constructGrid();
        }
    }

    /*
         public boolean deletePick()
         {
             if(patches == null) return false;
             StsPickPatch patch = (StsPickPatch)patches.getLast();
             if(patch == null) return false;
             patches.delete(patch);
             return true;
         }
      */

    public void display(StsGLPanel3d glPanel3d)
    {
        if (patches == null)
        {
            return;
        }
        patches.forEach("displaySeed", glPanel3d);
    }

    /*
         public void display(StsGLPanel3d glPanel3d)
         {
             if(wizard == null) return;
             if(wizard.getHorpick() != this) return;

             int nPickPatches = patches.getSize();
             patches.forEach("displaySeed", glPanel3d);

             GL gl = glPanel3d.getGL();

             if(displayState == DELETE)
             {
                 if(textureTiles != null)
                 {
                     textureTiles.deleteTextures(gl);
                     textureTiles.deleteSurface(gl);
                 }
                 displayState = NO_DRAW;
                 return;
             }
             if(textureTiles == null && !initializeSurface(gl)) return;
             if(displayState == NO_DRAW) return;

             if(displayState == BUILD)
             {
                 if(surface.getPointsZ() == null)
                 {
                     displayState = NO_DRAW;
                     return;
                 }
                 textureTiles.deleteSurface(gl);
                 textureTiles.constructSurface(surface, gl);
                 displayState = DRAW;
             }
             if(colorListChanged)
             {
                 gl.deleteLists(colorListNum, 1);
                 colorListNum = getColorListNum(gl);
                 colorListChanged = false;
             }

             gl.glDisable(GL.GL_LIGHTING);
             gl.glEnable(GL.GL_BLEND);
             gl.glEnable(GL.GL_TEXTURE_2D);

             gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
             gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
             gl.glShadeModel(GL.GL_FLAT);

             gl.glCallList(colorListNum);

             textureTiles.displayTiles(gl, true, patchIndexes, addData);
             addData = false;

             gl.glDisable(GL.GL_TEXTURE_2D);
             gl.glDisable(GL.GL_BLEND);
             gl.glEnable(GL.GL_LIGHTING);
         }
      */
    public void display2d(GL gl, int currentDirNo, float dirCoordinate, boolean axesFlipped)
    {
        if (patches == null)
        {
            return;
        }
        try
        {
            int nPatches = patches.getSize();
            gl.glDisable(GL.GL_LIGHTING);
            for (int n = 0; n < nPatches; n++)
            {
                StsPickPatch patch = (StsPickPatch) patches.getElement(n);
                StsColor stsColor = patch.getStsColor();
                StsGridPoint seedPoint = patch.getSeedPoint();
                if (seedPoint == null) continue;

                float[] xyz = seedPoint.getXYZorT();
                float x, y;
                switch (currentDirNo)
                {
                    case StsSeismicVolume.XDIR:
                        if (StsMath.sameAsTol(xyz[0], dirCoordinate, seismicVolume.getXInc()))
                        {
                            displaySeed2d(gl, xyz[1], xyz[2], stsColor);
                        }
                        break;
                    case StsSeismicVolume.YDIR:
                        if (StsMath.sameAsTol(xyz[1], dirCoordinate, seismicVolume.getYInc()))
                        {
                            displaySeed2d(gl, xyz[0], xyz[2], stsColor);
                        }
                        break;
                    case StsSeismicVolume.ZDIR:
                        // don't check because the seed point is adjusted vertically and may be off the zPlane
//                    if(StsMath.sameAsTol(xyz[2], dirCoordinate, line2d.getZInc()))
                    {
                        if (!axesFlipped)
                        {
                            displaySeed2d(gl, xyz[0], xyz[1], stsColor);
                        }
                        else
                        {
                            displaySeed2d(gl, xyz[1], xyz[0], stsColor);
                        }
                    }
                    break;
                }
            }
        }
        catch (Exception e)
        {
            StsException.systemError("StsHorpick.display2d() failed.");
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);

        }
    }

    public void displaySeed2d(GL gl, float x, float y, StsColor stsColor)
    {
//		currentModel.glPanel3d.setViewShift(gl, com.Sts.Utilities.StsGraphicParameters.vertexShift);
//		double viewShift = currentModel.glPanel3d.currentShift;
//		System.out.println("Current view shift: " + viewShift);
        StsGLDraw.drawPoint2d(x, y, StsColor.BLACK, gl, 8);
        StsGLDraw.drawPoint2d(x, y, stsColor, gl, 4);
//		currentModel.glPanel3d.resetViewShift(gl);
    }

    /*
         public boolean initializeSurface()
         {
             if(surface == null || patchIndexes == null) return false;
             surface.setNewTextureDisplayable(this);
             surface.setTextureChanged(true);
             surface.setSurfaceChanged(true);
             return true;
         }
      */
    /*
         public void clearTextureDisplay()
         {
             textureTiles.setDeleteTextures();
         }
      */
    // following required for StsTreeObjectI interface
    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsHorpick.class, "name", true, "Name:"),
                new StsColorComboBoxFieldBean(StsHorpick.class, "stsColor", "Color:", currentModel.getSpectrum("Basic").getStsColors())
            };
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if (propertyFields != null) return propertyFields;
        textureDisplayableListBean = new StsComboBoxFieldBean(StsHorpick.class, "surfaceTexture", "Property");
        colorscaleBean = new StsEditableColorscaleFieldBean(StsHorpick.class, "colorscale");
        propertyFields = new StsFieldBean[]
            {
                textureDisplayableListBean,
                new StsFloatFieldBean(StsHorpick.class, "minCorrelFilter", 0.0f, 1.0f, "Min Correl Filter"),
                new StsBooleanFieldBean(StsHorpick.class, "applyMinCorrelFilter", false, "Apply min correl filter"),
                colorscaleBean
            };
        return propertyFields;
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
        checkEnableCorrelCoefBean();
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.setCurrentObject(this);
//        currentModel.glPanel3d.checkAddView(StsViewXP.class);
//        currentModel.glPanel3d.cursor3d.setCurrentCrossplot(this);
        currentModel.win3dDisplayAll();
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsSeismicVolume getSeismicVolume()
    {
        return seismicVolume;
    }

    public void setSeismicVolume(StsSeismicVolume seismicVolume)
    {
        this.seismicVolume = seismicVolume;
    }

    public StsSurface getSurface()
    {
        return surface;
    }

    public void setSurface(StsSurface surface)
    {
        this.surface = surface;
        checkInitializeArrays();
        initializeSurfaceTextures();
        dbFieldChanged("surface", surface);
    }

    public void checkInitializeArrays()
    {
        if (patchIndexes != null && corCoefs != null)
        {
            return;
        }
        initializeArrays();
    }

    private void initializeArrays()
    {
        int nRows = surface.getNRows();
        int nCols = surface.getNCols();
        int nPoints = nRows * nCols;
        patchIndexes = new byte[nPoints];
        for (int n = 0; n < nPoints; n++)
        {
            patchIndexes[n] = -1;
        }
        corCoefs = new float[nRows][nCols];
    }

    private void initializeSurfaceTextures()
    {
        correlCoefsSurfaceTexture = new CorrelCoefsSurfaceTexture(this);
        patchColorSurfaceTexture = new PatchColorSurfaceTexture(this);
    }

    public void deleteArrays()
    {
        patchIndexes = null;
        corCoefs = null;
    }

    public void cleanup()
    {
        wizard = null;
    }

    public boolean isComplete()
    {
        return surface != null;
    }

    public boolean delete()
    {
        if (surface != null)
        {
            surface.delete();
        }
        super.delete();
        if (wizard != null)
        {
            wizard.deleteHorpick();
        }
        return true;
    }

    public void setIsIterative(boolean isIterative)
    {
        if (this.isIterative == isIterative) return;
        this.isIterative = isIterative;
        resetCorCoefColorscaleRange();
        wizard.runHorpick.panel.setIsIterative(isIterative);
    }

    private void resetCorCoefColorscaleRange()
    {
        float[] range = getCorCoefRange();
        if (!corCoefsColorscale.setRange(range)) return;
    }

    public boolean getIsIterative()
    {
        return isIterative;
    }

    public float getAutoCorInc()
    {
        return autoCorInc;
    }

    public float getAutoCorMax()
    {
        return autoCorMax;
    }

    public float getAutoCorMin()
    {
        return autoCorMin;
    }

    public void setAutoCorInc(float autoCorInc)
    {
        this.autoCorInc = autoCorInc;
    }

    public void setAutoCorMax(float autoCorMax)
    {
        this.autoCorMax = autoCorMax;
        resetCorCoefColorscaleRange();
    }

    public void setAutoCorMin(float autoCorMin)
    {
        this.autoCorMin = autoCorMin;
        resetCorCoefColorscaleRange();
    }

    public float getManualCorMin()
    {
        return manualCorMin;
    }

    public void setManualCorMin(float manualCorMin)
    {
        this.manualCorMin = manualCorMin;
        resetCorCoefColorscaleRange();
    }

    public byte[] getPatchIndexes()
    {
        return patchIndexes;
    }

    public int getPatchIndex(int row, int col)
    {
        int nCols = surface.getNCols();
        byte byteIndex = patchIndexes[row * nCols + col];
        return StsMath.signedByteToUnsignedInt(byteIndex);
    }

    public float[][] getCorCoefs()
    {
        return corCoefs;
    }

    public float getMinCorrel()
    {
        if (corCoefsColorscale == null)
        {
            return Math.min(autoCorMin, manualCorMin);
        }
        float[] editRange = corCoefsColorscale.getEditRange();
        return editRange[0];
    }

    public float getMaxCorrel()
    {
        if (corCoefsColorscale == null)
        {
            return 1.0f;
        }
        float[] editRange = corCoefsColorscale.getEditRange();
        return editRange[1];
    }

    public void setWizard(StsHorpickWizard wizard)
    {
        this.wizard = wizard;
    }

    public void printMemorySummary()
    {
        if (surface == null)
        {
            return;
        }
        int memoryUsed;
        int nPoints = surface.getNRows() * surface.getNCols();
        memoryUsed = 5 * nPoints;
        StsMessageFiles.infoMessage("Horpick arrays memory used: " + memoryUsed);
        memoryUsed = 24 * nPoints;
        StsMessageFiles.infoMessage("Surface displayList memory used: " + memoryUsed);
    }

    public void setDisplayPropertyName(String propertyName)
    {
        if (propertyName.equals(displayPropertyNone))
        {
            surface.setDefaultSurfaceTexture(null);
        }
        else if (propertyName.equals(displayPropertyPatchColor))
        {
            surface.setSurfaceTexture(patchColorSurfaceTexture, null);
        }
        else if (propertyName.equals(displayPropertyCorrelCoefs))
        {
            surface.setSurfaceTexture(correlCoefsSurfaceTexture, null);
        }
        else
        {
            StsSeismicVolumeClass volumeClass = (StsSeismicVolumeClass) currentModel.getCreateStsClass(StsSeismicVolume.class);
            StsSeismicVolume seismicVolume = (StsSeismicVolume) volumeClass.getObjectWithName(propertyName);
            surface.setSurfaceTexture(seismicVolume.getSurfaceTexture(surface), null);
        }
    }

    public StsSurfaceTexture getPatchColorSurfaceTexture()
    {
        return patchColorSurfaceTexture;
    }

    public StsSurfaceTexture getCorrelCoefsSurfaceTexture(StsSurface surface)
    {
        if (correlCoefsSurfaceTexture != null)
        {
            addCorCoefListeners(surface);
        }
        return correlCoefsSurfaceTexture;
    }

    public String getDisplayPropertyName()
    {
        return displayPropertyPatchColor;
    }

    public StsColorscale getCorCoefsColorscale()
    {
        if (corCoefsColorscale == null)
        {
            createCorCoefsColorscale();
        }
        return corCoefsColorscale;
    }

    public StsColorscalePanel getCorCoefsColorscalePanel()
    {
        return correlCoefsSurfaceTexture.getColorscalePanel();
    }

//    public void setCorCoefsColorscalePanel(StsColorscalePanel colorscalePanel) { correlCoefsTextureDisplayable.setColorscalePanel(colorscalePanel); }

    public StsSpectrum getSpectrum(String name)
    {
        return currentModel.getSpectrum(name);
    }

    public float getCorrelCoefMinFilter()
    {
        return minCorrelFilter;
    }

    public void setMinCorrelFilter(float minCorrelFilter)
    {
        if (this.minCorrelFilter == minCorrelFilter) return;
//		surface.applyMinCorrelFilter(applyMinCorrelFilter, correlCoefTextureDisplayable);
        StsModelSurfaceClass modelSurfaceClass = (StsModelSurfaceClass) currentModel.getCreateStsClass(StsModelSurface.class);
        StsModelSurface modelSurface = modelSurfaceClass.getModelSurfaceForOriginal(surface);
        if (modelSurface != null)
        {
            this.minCorrelFilter = minCorrelFilter;
            modelSurface.applyMinCorrelFilter(true, correlCoefsSurfaceTexture);
            currentModel.win3dDisplay();
        }
    }

    public float getMinCorrelFilter()
    {
        return Math.max(minCorrelFilter, getMinCorrel());
    }

    public StsColorscale getColorscale()
    {
        if (surfaceTexture == correlCoefsSurfaceTexture)
        {
            return getCorCoefsColorscale();
        }
        else
        {
            return null;
        }
    }

    private void createCorCoefsColorscale()
    {
        if(corCoefsColorscale != null) return;
        StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
        float[] range = getCorCoefRange();
        corCoefsColorscale = new StsColorscale("Correl Coefs", spectrum, range);
        addCorCoefListeners(surface);
    }

    public void addCorCoefListeners(StsSurface surface)
    {
        corCoefsColorscale.addActionListener(surface);
        corCoefsColorscale.addItemListener(surface);
    }

    private float[] getCorCoefRange()
    {
        if (isIterative)
            return new float[]{autoCorMin, 1.0f};
        else
            return new float[]{manualCorMin, 1.0f};
    }

    public void setColorscale(StsColorscale colorscale)
    {
        if (colorscale != corCoefsColorscale) return;
        if (surface == null) return;
        surface.setTextureChanged();
        surface.colorListChanged();
        currentModel.win3dDisplayAll();
    }

    public boolean getApplyMinCorrelFilter()
    {
        return applyMinCorrelFilter;
    }

    public void setApplyMinCorrelFilter(boolean applyMinCorrelFilter)
    {
        if (this.applyMinCorrelFilter == applyMinCorrelFilter) return;
        this.applyMinCorrelFilter = applyMinCorrelFilter;
        surface.setTextureChanged(correlCoefsSurfaceTexture);
        currentModel.win3dDisplay();
    }

    class PatchColorSurfaceTexture extends StsSurfaceTexture
    {
        StsHorpick horpick;
        int colorDisplayListNum = 0;

        PatchColorSurfaceTexture(StsHorpick horpick)
        {
            super(horpick.surface);
            this.horpick = horpick;
        }

        /** if horpick is deleted, it's index will now be -1; classes display this texture should delete it from their list */
        public boolean isDisplayable()
        {
            return horpick != null && horpick.isPersistent();
        }

        public byte[] getTextureData()
        {
            byte[] patchIndexes = horpick.getPatchIndexes();
            if (horpick.getMinCorrel() >= horpick.getMinCorrelFilter())
            {
                return patchIndexes;
            }

            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            byte[] textureData = new byte[nRows * nCols];
            System.arraycopy(patchIndexes, 0, textureData, 0, nRows * nCols);
            byte[][] surfacePointsNull = surface.getPointsNull();
            int n = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++, n++)
                {
                    if (surfacePointsNull[row][col] != StsSurface.SURF_PNT)
                    {
                        textureData[n] = -1;
                    }
                }
            }
            return textureData;
        }

        public String getName()
        {
            return StsHorpick.displayPropertyPatchColor;
        }

        public String toString()
        {
            return getName();
        }

        public StsColorscale getColorscale()
        {
            return null;
        }

        public float[] getHistogram()
        {
            return null;
        }

        public void deleteColorDisplayList(GL gl)
        {
            if (colorDisplayListNum == 0)
            {
                return;
            }
            gl.glDeleteLists(colorDisplayListNum, 1);
            colorDisplayListNum = 0;
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            deleteColorDisplayList(gl);

            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsHorpick.getColorListNum(): Failed to allocate a display list");
                return 0;
            }

            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
            arrayRGBA = null;
        }

        /*
                public void createColorTLUT(GL gl, boolean nullsFilled)
                {
                    float[][] arrayRGBA = computeRGBAArray(nullsFilled);
                    StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
                }
        */
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }

        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            //		StsSpectrum spectrum = horpick.getSpectrum("Basic");
            float[][] arrayRGBA = new float[4][256];
            StsObjectRefList patches = horpick.getPatches();
            int nPatches = patches.getSize();
            for (int n = 0; n < nPatches; n++)
            {
                StsPickPatch patch = (StsPickPatch) patches.getElement(n);
                StsColor color = patch.getStsColor();
                int index = patch.getIndex();
                if (index == -1) continue;  // shouldn't happen:  find bug!
                float[] rgba = color.getRGBA();
                for (int i = 0; i < 4; i++)
                {
                    arrayRGBA[i][index] = rgba[i];
                }
            }
            return arrayRGBA;
        }
    }

    class CorrelCoefsSurfaceTexture extends StsSurfaceTexture
    {
        StsHorpick horpick;
        int colorDisplayListNum = 0;
        StsColorscalePanel colorscalePanel;

        CorrelCoefsSurfaceTexture(StsHorpick horpick)
        {
            super(horpick.surface);
            this.horpick = horpick;
            createCorCoefsColorscale();
            addCorCoefListeners(surface);
        }

        /** if horpick is deleted, it's index will now be -1; classes display this texture should delete it from their list */
        public boolean isDisplayable()
        {
            return horpick != null && horpick.isPersistent();
        }

        public byte[] getTextureData()
        {
            float[][] corCoefs = horpick.getCorCoefs();
            StsSurface surface = horpick.surface;
            byte[][] surfacePointsNull = surface.getPointsNull();

            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            byte[] textureData = new byte[nRows * nCols];
            float[] range = getColorscale().getRange();
            float corMin = range[0];

            float scale = StsMath.floatToUnsignedByteScale(range);
            float scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, range[0]);

            int n = 0;
            for (int row = 0; row < nRows; row++)
            {
                for (int col = 0; col < nCols; col++)
                {
                    float corCoef = corCoefs[row][col];
                    if (surfacePointsNull[row][col] != StsSurface.SURF_PNT)
                        textureData[n++] = -1;
                    else if (horpick.applyMinCorrelFilter && horpick.corCoefs[row][col] < horpick.minCorrelFilter)
                        textureData[n++] = -1;
                    else
                        textureData[n++] = StsMath.floatToUnsignedByte254WithScale(corCoef, scale, scaleOffset);
                }
            }
            return textureData;
        }

        public float getDataMin()
        {
            if (isIterative)
                return autoCorMin;
            else
                return manualCorMin;
        }

        public float getDataMax()
        {
            return 1.0f;
        }

        public String getName()
        {
            return StsHorpick.displayPropertyCorrelCoefs;
        }

        public String toString()
        {
            return getName();
        }

        public StsColorscale getColorscale()
        {
            return horpick.getCorCoefsColorscale();
        }

        public float[] getHistogram()
        {
            if (colorscalePanel == null) getColorscalePanel();
            return colorscalePanel.getHistogram();
        }

        public void deleteColorDisplayList(GL gl)
        {
            if (colorDisplayListNum == 0)
            {
                return;
            }
            gl.glDeleteLists(colorDisplayListNum, 1);
            colorDisplayListNum = 0;
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            deleteColorDisplayList(gl);

            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsHorpick.getColorListNum(): Failed to allocate a display list");
                return 0;
            }

            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            //        if(mainDebug) System.out.println("Color 0: " + arrayRGBA[0][0] + " "  + arrayRGBA[1][0] + " "+ arrayRGBA[2][0] + " "+ arrayRGBA[3][0]);
            arrayRGBA = null;
        }

        /*
                public void createColorTLUT(GL gl, boolean nullsFilled)
                {
                    float[][] arrayRGBA = computeRGBAArray(nullsFilled);
                    StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
                }
        */
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }

        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            getColorscalePanel();
            Color[] colors = colorscalePanel.getColorscale().getNewColorsInclTransparency();
            if (nullsFilled)
                colors[colors.length - 1] = new Color(colors[colors.length - 1].getRed(), colors[colors.length - 1].getGreen(),
                    colors[colors.length - 1].getBlue(), 255);
            int nColors = colors.length;
            float[][] arrayRGBA = new float[4][nColors];
            float[] rgba = new float[4];
            for (int n = 0; n < nColors; n++)
            {
                colors[n].getComponents(rgba);
                for (int i = 0; i < 4; i++)
                {
                    arrayRGBA[i][n] = rgba[i];
                }
            }
            return arrayRGBA;
        }

        public StsColorscalePanel getColorscalePanel()
        {
            if (colorscalePanel != null) return colorscalePanel;

            StsColorscale colorscale = horpick.getCorCoefsColorscale();
            colorscalePanel = new StsColorscalePanel(colorscale, StsParameters.HORIZONTAL, false, false);
            colorscalePanel.setMinimumSize(new Dimension(200, 350));
            colorscalePanel.setPreferredSize(new Dimension(200, 350));
            colorscalePanel.addActionListener(horpick.getSurface());
            return colorscalePanel;
        }
    }
}