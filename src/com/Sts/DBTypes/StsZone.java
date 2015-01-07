//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsZone extends StsGridBoundingBox implements StsSelectable, StsTreeObjectI // , StsSurfaceTextureDisplayable
{
    protected StsSurface topSurface;
    protected StsSurface botSurface;
    protected StsModelSurface topModelSurface;
    protected StsModelSurface baseModelSurface;
    protected StsObjectRefList subZones;

    protected StsColor stsColor;
    protected boolean displayGrid = true;
    protected boolean displayFill = true;
    protected boolean displayGaps = true;

    transient protected StsList zoneBlocks;
    transient protected StsList subZoneSurfaces;
    transient private ArrayList<StsPropertyVolume> zoneProperties = new ArrayList<StsPropertyVolume>();

    /** Zone consists of a set of subZones. Each subZone has a specified number of layers.
     *  The number of layerGrids is one more than the sum all the layers in all the subZones.
     */
    /** sequential zone number assigned in building the model; numbered from to down beginning with zero */
    transient protected int zoneNumber;
    /** beginning layerNo at top of zone; layers are numbered from zero at the top of the model sequentially down through each of the zones */
    transient protected int topLayerNumber;
    /**
     * Layers and layer grids are numbered in sequence from the top of the model down. If there are a total of N layers,
     * there will be N+1 layerGrids.  The layerGrid above layer n is n, and the layerGrid below it is n+1.
     * SubZoneBotLayerGridNumbers are the layerGrid numbers at the bottom of each subZone, i.e., the number of the layer grid
     * just below the bottom layer in the subZone.
     */
    transient protected int[] subZoneBotLayerGridNumbers;

    protected int status = UNDEFINED;
    // protected int nSubZones = DEFAULT_N_SUBZONES; // obsolete
    protected int subZoneType = SUBZONE_UNIFORM;
    protected float subZoneThickness;
    protected StsWellZoneSet wellZoneSet = null;

    protected boolean hasSubHorizonSurfaces = false;

    // zone base and top surfaces and subHorizon surfaces in between ordered bot to top
    protected StsObjectRefList subHorizonSurfaces = null;
/*
    protected StsPropertyVolumeOld propertyVolume;
    protected int propVolBotLayer;
    protected int propVolTopLayer;
*/
    /** Source of texture for surface currently being used. */
    //    transient StsSurfaceTextureDisplayable newTextureDisplayable = this;
    /** Texture has been changed: replace texture subImage2D with new texture */
    transient boolean textureChanged = false;

    transient protected float[][] thickness = null;
    transient protected byte[][] nullTypes = null;

    static protected int currentLayer = 1;
    static StsZoneClass zoneClass = null;
    static protected StsObjectPanel objectPanel = null;

    // *** comment out types not yet implemented! ***
    static public final String[] subZoneTypeStrings =
        {"Stratigraphic", "Onlap", "Offlap"
        };

    static StsComboBoxFieldBean textureDisplayableListBean;
    static public StsFieldBean[] displayFields = null;

    static StsEditableColorscaleFieldBean colorscaleBean;
    static public StsFieldBean[] propertyFields = null;

    // constants
    static public final int SUBZONE_UNIFORM = 0;
    static public final int SUBZONE_ONLAP = 1;
    static public final int SUBZONE_OFFLAP = 2;
    static public final int SUBZONE_TRUNCATED = 3;

    static public final int SUBZONE_MAX = 500;

    static public final int DEFAULT_N_SUBZONES = 10;

    static public final String ZONE_SURFACE_PREFIX = "Zone-";
    static public final String ZONE_SURFACE_TOP = "-top";
    static public final String ZONE_SURFACE_BASE = "-base";

    static public final int UNDEFINED = StsParameters.UNDEFINED;
    static public final int NOT_BUILDABLE = 1;
    static public final int BUILDABLE = 2;
    static public final int BUILT = 3;

    // Convenience copies of flags
    static final float largeFloat = StsParameters.largeFloat;
    static final float nullZValue = StsParameters.nullValue;
    static public final byte SURF_PNT = StsGridPoint.SURF_PNT;
    static public final byte SURF_BOUNDARY = StsGridPoint.SURF_BOUNDARY;
    static public final byte SURF_GAP = StsGridPoint.SURF_GAP;
    static public final byte SURF_GAP_SET = StsGridPoint.SURF_GAP_SET;

    static final boolean debug = true;
    
    /** constructor  for surfaces */
    public StsZone(String name, StsModelSurface topSurface, StsModelSurface botSurface) throws StsException
    {
        super(false);
        setTopModelSurface(topSurface);
        setBaseModelSurface(botSurface);
        addBoundingBox(topSurface);
        addBoundingBox(botSurface);
        if (name == null)
        {
            name = topSurface.getName();
        }
        setName(name);

        stsColor = topSurface.getStsColor();
//        constructHardWiredSubZones();
        // constructSubZones();
        // constructProperties();
        addToModel();
//		refreshObjectPanel();
    }

    /** constructor for well zones */
    public StsZone(StsWellZoneSet wellZoneSet) throws StsException
    {
        if (wellZoneSet == null)
        {
            throw new StsException(StsException.FATAL, "StsZone.setTop:" + " Cannot create a zone with a null well zone set");
        }
        setWellZoneSet(wellZoneSet);
        topSurface = null;
        botSurface = null;
        setName(wellZoneSet.getName());
        stsColor = wellZoneSet.getStsColor();
        wellZoneSet.setParentZone(this);
    }

    /** DB constructor */
    public StsZone()
    {
    }

    public void constructHardWiredSubZones(int nSubZones)
    {
        constructSubZones(nSubZones);
    }

    public String getName()
    {
        return super.getName();
    }

    public boolean initialize(StsModel model)
    {
        if (zoneClass == null)
        {
            zoneClass = (StsZoneClass) model.getCreateStsClass(this);
        }
        computeNLayers();
        return true;
    }

    public void clearTransientArrays()
    {
        zoneProperties = new ArrayList<StsPropertyVolume>();
    }

    /*
      public void initTextureDisplayableList()
      {
       if(textureDisplayableList != null) return;

       StsSeismicVolume[] seismicVolumes = (StsSeismicVolume[])currentModel.getCastObjectList(StsSeismicVolume.class);
       int nListItems = seismicVolumes.length;

       textureDisplayableList = new StsSurfaceTextureDisplayable[nListItems];
       for(int n = 0; n < seismicVolumes.length; n++)
        textureDisplayableList[n] = seismicVolumes[n];

       textureDisplayableListBean.setListItems(textureDisplayableList);
      }

      public StsSurfaceTextureDisplayable getNewTextureDisplayable()
      {
       initTextureDisplayableList();
       return newTextureDisplayable;
      }
      */
    public void display(StsGLPanel3d glPanel3d, String displayModeString)
    {
        if (!isVisible) return;
        if (zoneBlocks == null) return;

        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            StsBlock block = zoneBlock.getBlock();
            if (!block.isVisible) continue;
            StsPropertyVolume propertyVolume = block.getCurrentPropertyVolume();
            zoneBlock.display(currentModel, glPanel3d, displayModeString, propertyVolume);
        }
    }

    public void displayLayer(StsGLPanel3d glPanel3d, int layer, boolean displayGaps)
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            zoneBlock.displayGridSlice(glPanel3d, layer, displayGaps);
        }
    }

    // Accessors
    public void setTopSurface(StsModelSurface topSurface) throws StsException
    {
        if (topSurface == null)
        {
            throw new StsException(StsException.FATAL, "StsZone.setTop:" + " Cannot create a zone with a null top surface");
        }

        if (this.topSurface != topSurface)
        {
            this.topSurface = topSurface;
            topSurface.setZoneBelow(this);

//            if(botSurface != null)
//                setThicknessAndNullTypes(currentModel);
        }
    }

    public StsSurface getTopSurface()
    {
        return topSurface;
    }

    public void setBotSurface(StsModelSurface botSurface) throws StsException
    {
        if (botSurface == null)
        {
            throw new StsException(StsException.FATAL, "StsZone.setBotSurface:" + " Cannot create a zone with a null base surface");
        }
        if (this.botSurface != botSurface)
        {
            this.botSurface = botSurface;
            botSurface.setZoneAbove(this);

//            if(topSurface != null)
//                setThicknessAndNullTypes(currentModel);
        }
    }

    public StsSurface getBotSurface()
    {
        return botSurface;
    }

    public boolean subZoneChanged(StsSubZone subZone)
    {
        computeNLayers();
        if (zoneBlocks != null)
        {
            zoneBlocks.forEach("reconstructSubZoneGrid", subZone);
        }
        /*
        if(propertyVolume != null && propertyVolume.getNLayers() != nLayerGrids)
        {
            propertyVolume = null;
        }
        */
        deleteDisplayLists();
        return true;
    }

    public void displaySubZoneColorsChanged()
    {
        deleteDisplayLists();
    }

    private void computeNLayers()
    {
        if (subZones == null) return;
        int nSubZones = subZones.getSize();
        int nLayers = 0;
        subZoneBotLayerGridNumbers = new int[nSubZones];

        for (int n = 0; n < nSubZones; n++)
        {
            StsSubZone subZone = (StsSubZone) subZones.getElement(n);
            nLayers += subZone.getNLayers();
            subZoneBotLayerGridNumbers[n] = nLayers;
        }
        setNLayers(nLayers);
    }

    // TODO doesn't look right, but not currently active; so check carefully if you need it
    public int[] getSubZoneAndLayer(int zoneLayer)
    {
        if (zoneLayer < 0 || zoneLayer >= sliceMax) return null;

        int nSubZones = subZones.getSize();
        if (nSubZones == 0) return null;

        int nBotLayer = -1, nTopLayer = 0;
        for (int n = 0; n < nSubZones; n++)
        {
            StsSubZone subZone = (StsSubZone) subZones.getElement(n);
            nTopLayer = nBotLayer + 1;
            nBotLayer += subZone.nLayers;
            if (zoneLayer <= nBotLayer)
            {
                return new int[]{n, zoneLayer - nTopLayer};
            }
        }
        return null;
    }

    public void setNLayers(int nLayers)
    {
        sliceMin = topLayerNumber;
        sliceMax = sliceMin + nLayers;
    }

    public int checkLimitLayerNumber(int nLayer)
    {
        int botLayerNumber = topLayerNumber + sliceMax - sliceMin;
        if(debug && (nLayer  < topLayerNumber || nLayer > botLayerNumber))
            StsException.systemDebug(this, "checkLimitLayerNumber", "nLayerF: " + nLayer + " is out of range " + topLayerNumber + " - " + botLayerNumber);

        if(nLayer < topLayerNumber)
            return topLayerNumber;
        else if(nLayer > botLayerNumber)
            return botLayerNumber;
        else
            return nLayer;
    }


    public int checkLimitLayerNumber(float nLayerF)
    {
        int botLayerNumber = topLayerNumber + sliceMax - sliceMin;
        if(debug && (nLayerF  < topLayerNumber || nLayerF > botLayerNumber))
            StsException.systemDebug(this, "checkLimitLayerNumber", "nLayerF: " + nLayerF + " is out of range " + topLayerNumber + " - " + botLayerNumber);

        if(nLayerF < topLayerNumber)
            return topLayerNumber;
        else if(nLayerF > botLayerNumber)
            return botLayerNumber;
        else
            return (int)nLayerF;
    }

    public void setNSubZones(int nSubZones)
    {
        if (subZones == null || subZones.getSize() == nSubZones) return;
        nSubZonesChanged(nSubZones);
    }

    /*
      static public BoundedRangeModel getLayerRange()
      {
       int nLayers = getTotalNLayers();
       if(nLayers == 0) return null;
       DefaultBoundedRangeModel rangeModel = new DefaultBoundedRangeModel();
       try
       {
        rangeModel.setRangeProperties(currentLayer, 0, 0, nLayers, false);
       }
       catch(Exception e) { e.printStackTrace(); }
       return rangeModel;
      }

      static public int getTotalNLayers()
      {
       try
       {
        StsClass zones = currentModel.getCreateStsClass(StsZone.class);
        int nZones = zones.getSize();
        int nLayers = 0;
        for (int n = 0; n < nZones; n++)
        {
         StsZone zone = (StsZone)zones.getElement(n);
         nLayers += zone.getNLayers();
        }
        return nLayers;
       }
       catch (Exception e) { return 0; }
      }

      static public void setLayerRange(BoundedRangeModel range)
      {
       currentLayer = range.getValue();
      }

      static public int getCurrentLayer() { return currentLayer; }
      static public void setCurrentLayer(int layer) { currentLayer = layer; }
      */
    private void nSubZonesChanged(int nSubZones)
    {
//        if (wellZoneSet != null && !wellZoneSet.setNSubZones(nSubZones, false)) return false;

//        if(propertyVolume != null && propertyVolume.getNSubZones() != nSubZones)
//        propertyVolume = null;

        deleteDisplayLists();
        constructSubZones(nSubZones);
        if (zoneBlocks != null)
        {
            zoneBlocks.forEach("constructGrid", null);
        }
    }

    /*  nSubZones now fixed by number of intermediate horizons
      public boolean setNSubZones(int nSubZones)
      {
       if (nSubZones < 1 || nSubZones > SUBZONE_MAX) return false;
       this.nSubZones = nSubZones;
       if (wellZoneSet != null && !wellZoneSet.setNSubZones(nSubZones, false)) return false;
       if(zoneBlocks != null) zoneBlocks.forEach("reconstruct3dGrid", null);
       if(propertyVolume != null && propertyVolume.getNSubZones() != nSubZones)
        propertyVolume = null;
       deleteDisplayLists();
       return true;
      }
      */
    public int getNLayerGrids()
    {
        return sliceMax - sliceMin + 1;
    }

    public int getNLayers()
    {
        return sliceMax - sliceMin;
    }

    public StsObjectRefList getSubZones()
    {
        return subZones;
    }

    public int getNSubZones()
    {
        return subZones == null ? 0 : subZones.getSize();
    }

    public float getLayerF(float fraction)
    {
        int nLayers = getNLayers();
        return topLayerNumber + fraction * nLayers;
    }

    public boolean hasLayer(int nLayer)
    {
        return sliceMin <= nLayer && sliceMax - 1 >= nLayer;
    }

    public boolean hasLayerGrid(int nLayerGrid)
    {

        return sliceMin <= nLayerGrid && sliceMax >= nLayerGrid;
    }

    public StsColor getLayerColor(int layer)
    {
        return currentModel.getSpectrum("Basic").getColor(layer);
    }

    // comment out types not yet implemented!
    public boolean setSubZoneType(int type)
    {
        switch (type)
        {
            case SUBZONE_UNIFORM:

                /*
                        case SUBZONE_TRUNCATED:
                        case SUBZONE_ONLAP_OFFLAP:
                     */
                subZoneType = type;
                return true;
        }
        return false;
    }

    public String getSubZoneTypeStr()
    {
        return subZoneTypeStrings[subZoneType];
    }

    public boolean setSubZoneTypeStr(String type)
    {
        for (int i = 0; i < subZoneTypeStrings.length; i++)
        {
            if (type.equals(subZoneTypeStrings[i]))
            {
                subZoneType = i;
                return true;
            }
        }
        return false;
    }

    public int getSubZoneType()
    {
        return subZoneType;
    }

    public void setStsColor(StsColor color)
    {
        if (stsColor == color)
        {
            return; // no change
        }
        stsColor = color;
        if (wellZoneSet != null)
        {
            wellZoneSet.setStsColor(stsColor, true);
        }
        if (topModelSurface != null)
        {
            topModelSurface.setStsColor(color);
        }
        if (baseModelSurface != null && baseModelSurface.getZoneBelow() == null)
        {
            baseModelSurface.setStsColor(color);
        }
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setIsVisible(boolean b)
    {
        isVisible = b;
        if (b && currentModel != null)
        {
            zoneClass.setDisplayZones(true);
        }
    }

    public boolean getIsVisible()
    {
        return isVisible;
    }

    public boolean getDisplayGrid()
    {
        return displayGrid;
    }

    public boolean getDisplayGaps()
    {
        return displayGaps;
    }

    public boolean getDisplayFill()
    {
        return displayFill;
    }

    public void setWellZoneSet(StsWellZoneSet set)
    {
        wellZoneSet = set;
    }

    public StsWellZoneSet getWellZoneSet()
    {
        return wellZoneSet;
    }

    public StsModelSurface getTopModelSurface()
    {
        return topModelSurface;
    }

    public void setTopModelSurface(StsModelSurface surface)
    {
        topModelSurface = surface;
        topModelSurface.setZoneBelow(this);
    }

    public StsModelSurface getBaseModelSurface()
    {
        return baseModelSurface;
    }

    public void setBaseModelSurface(StsModelSurface surface)
    {
        baseModelSurface = surface;
        baseModelSurface.setZoneAbove(this);
    }

    /*
    public StsPropertyVolumeOld getPropertyVolume()
    {
        return propertyVolume;
    }
    */
    public void setSubZoneThickness(float subZoneThickness)
    {
        this.subZoneThickness = subZoneThickness;
    }

    public float getSubZoneThickness()
    {
        return subZoneThickness;
    }

    public StsObjectRefList getSubHorizonSurfaces()
    {
        return subHorizonSurfaces;
    }

    public StsZoneBlock[] getZoneBlocks()
    {
        return (StsZoneBlock[])zoneBlocks.getCastList(StsZoneBlock.class);
    }

    public void setDisplayFill(boolean b)
    {
        displayFill = b;
//        if(topModelSurface != null) topModelSurface.setDisplayFill(b);
//        if(baseModelSurface != null) baseModelSurface.setDisplayFill(b);
    }

    public void setDisplayGrid(boolean b)
    {
        displayGrid = b;
//        if(topModelSurface != null) topModelSurface.setDisplayGrid(b);
//        if(baseModelSurface != null) baseModelSurface.setDisplayGrid(b);
    }

    public void setDisplayGaps(boolean b)
    {
        displayGaps = b;
//        if(topModelSurface != null) topModelSurface.setDisplayGaps(b);
//        if(baseModelSurface != null) baseModelSurface.setDisplayGaps(b);
    }

    public void setStatus(int status)
    {
        this.status = status;
    }

    public int getStatus()
    {
        if (status == BUILT)
        {
            return status;
        }

        else if (topModelSurface != null && baseModelSurface != null)
        {
            status = BUILDABLE;
        }
        else if (topModelSurface != null && wellZoneSet != null)
        {
            status = BUILDABLE;
        }
        else if (baseModelSurface != null && wellZoneSet != null)
        {
            status = BUILDABLE;
        }
        else
        {
            status = NOT_BUILDABLE;

        }
        return status;
    }

    public void setThicknessAndNullTypes(StsModel model)
    {
        float[][] basePointsZ, topPointsZ;
        byte[][] basePointsNull, topPointsNull;
        float baseZ, topZ;

        basePointsZ = baseModelSurface.getPointsZ();
        topPointsZ = topModelSurface.getPointsZ();
        basePointsNull = baseModelSurface.getPointsNull();
        topPointsNull = topModelSurface.getPointsNull();

        int nRows = baseModelSurface.getNRows();
        int nCols = baseModelSurface.getNCols();
        int nPnts = nRows * nCols;

        thickness = new float[nRows][nCols];
        nullTypes = new byte[nRows][nCols];

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                topZ = topPointsZ[row][col];
                baseZ = basePointsZ[row][col];

                if (topZ > baseZ)
                {
                    topPointsZ[row][col] = baseZ;
                    thickness[row][col] = 0.0f;
                }
                else
                {
                    thickness[row][col] = baseZ - topZ;

                }
                byte topNullType = topPointsNull[row][col];
                byte baseNullType = basePointsNull[row][col];
                if (topNullType == SURF_PNT && baseNullType == SURF_PNT)
                {
                    nullTypes[row][col] = SURF_PNT;
                }
                else
                {
                    if (topNullType == SURF_BOUNDARY || baseNullType == SURF_BOUNDARY)
                    {
                        nullTypes[row][col] = SURF_BOUNDARY;
                    }
                    else if (topNullType == SURF_GAP || baseNullType == SURF_GAP)
                    {
                        nullTypes[row][col] = SURF_GAP;
                    }
                    else if (topNullType == SURF_GAP_SET || baseNullType == SURF_GAP_SET)
                    {
                        nullTypes[row][col] = SURF_GAP_SET;
                    }
                }
            }
        }
    }

    public StsList createModelSubZoneSurfaces(StsModel model)
    {
        int nSubZones = subZones.getSize();
        int nSubZoneSurfaces = nSubZones - 1;
        subZoneSurfaces = new StsList(nSubZoneSurfaces, 2);
        StsColor color = topModelSurface.getStsColor();

        float df = 1.0f / nSubZones;
        float f = df;
        for (int n = 0; n < nSubZoneSurfaces; n++)
        {
            StsModelSurface subZoneSurface = topModelSurface.cloneModelSurface();
            String name = new String("subZone-" + getName() + "-" + n);
            subZoneSurface.setName(name);
            subZoneSurface.setSubZoneNumber(n);
            subZoneSurface.setZone(this);
            subZoneSurfaces.add(subZoneSurface);
            interpolateSubZoneGrid(subZoneSurface, f);
            f += df; //    subZoneGrid.saveValuesToDataStore();
        }

        status = BUILT;

        return subZoneSurfaces;
    }

    private void interpolateSubZoneGrid(StsModelSurface subZoneSurface, float f)
    {
        int i, j, n;

        float[][] pointsZ, basePointsZ;
        byte[][] pointsNull, basePointsNull;

        pointsZ = subZoneSurface.getPointsZ();
        pointsNull = subZoneSurface.getPointsNull();

        basePointsZ = baseModelSurface.getPointsZ();
        basePointsNull = baseModelSurface.getPointsNull();

        int nRows = baseModelSurface.getNRows();
        int nCols = baseModelSurface.getNCols();
        int nPnts = nRows * nCols;

        for (int row = 0; row < nRows; row++)
        {
            for (int col = 0; col < nCols; col++)
            {
                pointsZ[row][col] = (basePointsZ[row][col] - f * thickness[row][col]);
                pointsNull[row][col] = nullTypes[row][col];
            }
        }
        subZoneSurface.constructGrid();
//		subZoneSurface.persistPoints();
    }

    /*
       public void setProperty(StsPropertyVolumeOld propVol, int propVolBotLayer, int propVolTopLayer)
       {
           if(baseModelSurface.getZoneBelow() == null)
           {
               baseModelSurface.setProperty(propVol, propVolBotLayer);
           }
           topModelSurface.setProperty(propVol, propVolTopLayer);

           propertyVolume = propVol;
           this.propVolBotLayer = propVolBotLayer;
           this.propVolTopLayer = propVolTopLayer;

           if(propVol != null)
           {
   //            setDisplayFill(true);
   //            setIsVisible(true);
           }
       }
    */
    public void deleteDisplayLists()
    {
        if (zoneBlocks != null)
        {
            zoneBlocks.forEach("deleteDisplayLists", null);
        }
    }

    /*
       public void setProperty(StsPropertyVolumeOld propVol)
       {
      int nPropSubZones = propVol.getNSubZones();

      baseModelSurface.setProperty(propVol, 0);
 //        propVol.setValuesOnSurface(baseModelSurface, 0);

      topModelSurface.setProperty(propVol, nPropSubZones-1);
 //        propVol.setValuesOnSurface(topModelSurface, nPropSubZones-1);

      if(subZoneSurfaces == null) return;

      int nSubZoneSurfaces = subZoneSurfaces.getSize();

      for(int n = 0; n < nSubZoneSurfaces; n++)
      {
       StsSurface subZoneSurface = (StsSurface)subZoneSurfaces.getElement(n);
       int nPropSubZone = Math.min(nPropSubZones-1, nPropSubZones*(n+1)/nSubZones);
       subZoneSurface.setProperty(propVol, nPropSubZone);
 //            propVol.setValuesOnSurface(subZoneSurface, nPropSubZone);
      }
       }
      */

    public boolean checkZoneSurfaces()
    {
        /*
           boolean baseChecked, topChecked;

           baseChecked = baseModelSurface != null &&
                baseModelSurface.isType(StsModelSurface.CHECKED);
           topChecked = topModelSurface != null &&
                topModelSurface.isType(StsModelSurface.CHECKED);

           if(baseChecked && topChecked) return true;

           if(topChecked)
           {
            if(baseModelSurface == null)
            {
             if(wellZoneSet != null)
             {
              baseModelSurface = getBaseFromTopAndWells(topModelSurface, wellZoneSet);
              return baseModelSurface != null;
             }
             else
             {
              baseModelSurface = getBaseFromTop(topModelSurface);
              return baseModelSurface != null;
             }
            }
            else
             return checkBaseSurface(topModelSurface);
           }
           if(baseChecked)
           {
            if(topModelSurface == null)
            {
             if(wellZoneSet != null)
             {
              topModelSurface = getTopFromBaseAndWells(baseModelSurface, wellZoneSet);
              return topModelSurface != null;
             }
             else
             {
              topModelSurface = getTopFromBase(baseModelSurface);
              return topModelSurface != null;
             }
            }
            else
             return checkTopSurface(topModelSurface);
           }
           */
        return true;
    }

    // for now construct proportional subZones
    public void constructSubZones(int nSubZones)
    {
        int n;
        float f, df;

        if (nSubZones < 1)
        {
            nSubZones = 1;
        }
        df = 1.0f / nSubZones;
        f = 0.0f;
        subZones = StsObjectRefList.constructor(nSubZones, 1, "subZones", this);
        for (n = 0; n < nSubZones; n++)
        {
            StsSubZone subZone = new StsSubZone(StsSubZone.SUBZONE_UNIFORM, f, f + df, this, n, true);
            subZones.add(subZone);
            f += df;
        }
        computeNLayers();
    }

    public StsPropertyVolume getPropertyVolume(StsPropertyType propertyType)
    {
        String propertyName = propertyType.name;
        for (StsPropertyVolume propertyVolume : zoneProperties)
        {
            if (propertyVolume.getName().equals(propertyName))
                return propertyVolume;
        }
        return null;
    }

    public StsPropertyVolume getPropertyVolume(String propertyName)
    {
        for (StsPropertyVolume propertyVolume : zoneProperties)
        {
            if (propertyVolume.propertyType.name.equals(propertyName))
                return propertyVolume;
        }
        return null;
    }

    public StsPropertyVolume getEclipsePropertyVolume(String eclipsePropertyName)
    {
        for (StsPropertyVolume propertyVolume : zoneProperties)
        {
            if (propertyVolume.propertyType.eclipseName.equals(eclipsePropertyName))
                return propertyVolume;
        }
        return null;
    }

    /*
    public void constructSubZones()
    {
     int n;

     StsSurface[] orderedSurfaces = currentModel.getOrderedImportedSurfaces();
     int nOrderedSurfaces = orderedSurfaces.length;

     int topSurfaceIndex = -1;
     int botSurfaceIndex = -1;
     for(n = 0; n < nOrderedSurfaces; n++)
     {
      if(topSurface == orderedSurfaces[n]) topSurfaceIndex = n;
      else if(botSurface == orderedSurfaces[n]) botSurfaceIndex = n;
     }

     // if we didn't find these surfaces (shouldn't happen) or there are no
     // subHorizons in between (bot and top are adjacent) then return
     if(topSurfaceIndex == -1 || botSurfaceIndex == -1) return;
     if(botSurfaceIndex - topSurfaceIndex <= 0) return;

     nSubZones = botSurfaceIndex - topSurfaceIndex;
     subHorizonSurfaces = new StsObjectRefList(nSubZones+1, 2);

     for(n = topSurfaceIndex; n <= botSurfaceIndex; n++)
     {
      StsSurface surface = orderedSurfaces[n];
      subHorizonSurfaces.add(surface);
     }

     subZones = new StsObjectRefList(nSubZones);
     StsSurface botSubZoneSurface, topSubZoneSurface;
     botSubZoneSurface = (StsSurface)subHorizonSurfaces.getElement(0);
     for(n = 0; n < nSubZones; n++)
     {
      topSubZoneSurface = botSubZoneSurface;
      botSubZoneSurface = (StsSurface)subHorizonSurfaces.getElement(n+1);
      StsSubZone subZone = new StsSubZone(botSubZoneSurface, topSubZoneSurface, this, n);
      subZones.add(subZone);
     }

     computeNLayers();
    }
    */
    public void checkLoadSubHorizonGrids()
    {
        int n;

        if (subHorizonSurfaces == null)
        {
            return;
        }

        int nSubHorizonSurfaces = subHorizonSurfaces.getSize();
        for (n = 0; n < nSubHorizonSurfaces; n++)
        {
            StsModelSurface surface = (StsModelSurface) subHorizonSurfaces.getElement(n);
            if (!surface.checkIsLoaded())
            {
                subHorizonSurfaces.delete();
                return;
            }
        }
    }

    public boolean buildZoneBlocks(StsBlock[] blocks, StsStatusUI status, float progress, float incProgress)
    {
        StsBlockGrid topBlockGrid, botBlockGrid;
        StsZoneBlock zoneBlock;

        try
        {
            if (blocks == null) return false;

            int nBlocks = blocks.length;
            incProgress = incProgress / nBlocks;
            zoneBlocks = new StsList(nBlocks, 1);

            for (StsBlock block : blocks)
            {
                topBlockGrid = block.getBlockGrid(topModelSurface);
                if (topBlockGrid == null)
                {
                    continue;
                }
                botBlockGrid = block.getBlockGrid(baseModelSurface);
                if (botBlockGrid == null)
                {
                    continue;
                }
                zoneBlock = new StsZoneBlock(this, block, botBlockGrid, topBlockGrid);
                zoneBlocks.add(zoneBlock);
                block.addBlockZone(zoneBlock);
                status.setProgress(progress += incProgress);
//                StsStatusArea.staticSetProgress(progress += incProgress);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsZone.buildZoneBlocks() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dir, float dirCoordinate, StsPoint startPoint, StsPoint endPoint, boolean isSliderDragging, String displayMode)
    {
        if (zoneBlocks == null)
        {
            return;
        }
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            zoneBlock.drawCursorSection(currentModel, glPanel3d, dir, dirCoordinate, startPoint, endPoint, isSliderDragging, displayMode);
        }
    }

    public boolean buildWellZoneSet(StsModel model)
    {
        StsModelSurface topSurface = getTopModelSurface();
        StsModelSurface botSurface = getBaseModelSurface();
        if (topSurface == null || botSurface == null)
        {
            return false;
        }

        StsMarker topMarker = topSurface.getMarker();
        StsMarker baseMarker = botSurface.getMarker();
        if (topMarker == null || baseMarker == null)
        {
            return false;
        }

        if (StsWellZoneSet.isNewWellZoneSet(model, topMarker, baseMarker))
        {
            StsWellZoneSet wzs = null;
            try
            {
                wzs = new StsWellZoneSet(StsWellZoneSet.STRAT, topMarker, baseMarker);
            }
            catch (StsException e)
            {
                return false;
            }
            wzs.setParentZone(this);
            setWellZoneSet(wzs);
            return true;
        }

        return false;
    }

    public StsZoneBlock getZoneBlockBelow(StsBlockGrid blockGrid)
    {
        if (zoneBlocks == null) return null;

        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (zoneBlock.getTopGrid() == blockGrid)
            {
                return zoneBlock;
            }
        }
        return null;
    }

    public StsZoneBlock getZoneBlockAbove(StsBlockGrid blockGrid)
    {
        if (zoneBlocks == null) return null;

        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (zoneBlock.getBottomGrid() == blockGrid)
            {
                return zoneBlock;
            }
        }
        return null;
    }

    /** delete any model-related objects */
     public boolean deleteModel()
     {
         zoneBlocks = null;
         return true;
     }

    public StsFieldBean[] getDisplayFields()
    {
        if (displayFields == null)
        {
            // textureDisplayableListBean = new StsComboBoxFieldBean(StsSurface.class, "newTextureDisplayable", "Property");
            displayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsZone.class, "isVisible", "Enable"),
                    new StsBooleanFieldBean(StsZone.class, "displayFill", "Fill"),
                    new StsBooleanFieldBean(StsZone.class, "displayGrid", "Grid")
                };
        }
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        /*
        if(propertyFields == null)
        {
            colorscaleBean = new StsEditableColorscaleFieldBean(StsSurface.class, "colorscale");
            propertyFields = new StsFieldBean[]
                    {
                            colorscaleBean
                    };
        }
        */
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
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        StsClass stsClass = currentModel.getCreateStsClass(StsZone.class);
        stsClass.selected(this);
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsColorscale getColorscale()
    {
        return null;
    }

    public int getZoneNumber()
    {
        return zoneNumber;
    }

    public void setZoneNumber(int nZone)
    {
        this.zoneNumber = nZone;
    }

    public int getTopLayerNumber()
    {
        return sliceMin;
    }

    public int getBottomLayerNumber()
    {
        return sliceMax - 1;
    }

    public void setTopLayerNumber(int topLayerNumber)
    {
        this.topLayerNumber = topLayerNumber;
        int nLayers = sliceMax - sliceMin;
        sliceMin = topLayerNumber;
        sliceMax = topLayerNumber + nLayers;
    }

    public boolean matches(String name, StsModelSurface topModelSurface, StsModelSurface baseModelSurface)
    {
        return this.topModelSurface == topModelSurface && this.baseModelSurface == baseModelSurface;

    }

    static public boolean constructLayerGrids(StsZone[] zones)
    {
        int nZones = zones.length;
        for (int n = 0; n < nZones; n++)
            if(!zones[n].constructLayerGrids()) return false;
        return true;
    }

    public boolean constructLayerGrids()
    {
        int nZoneBlocks = zoneBlocks.getSize();
        for (int n = 0; n < nZoneBlocks; n++)
        {
            StsZoneBlock zoneBlock = (StsZoneBlock) zoneBlocks.getElement(n);
            if (!zoneBlock.constructLayerGrids()) return false;
        }
        return true;
    }

    public ArrayList<StsPropertyVolume> getZoneProperties()
    {
        return zoneProperties;
    }
}