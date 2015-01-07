package com.Sts.MVC;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: StsProject is the class that maintains all the boundary conditions for the Project.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.DataCube.*;
import com.Sts.Utilities.DateTime.*;
import com.Sts.Utilities.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.text.*;
import java.util.*;


public class StsProject extends StsObject implements ActionListener, Serializable, StsTreeObjectI, StsCultureDisplayable
{
    public static final float AUTO = -1.0f;
    /** rotated box bounding surfaces and seismic which defines the 3 cursor planes */
    public StsRotatedGridBoundingBox rotatedBoundingBox;
    /** Unrotated box tightly bounding all project objects. This boundingBox carries the common origin which is shared with the rotatedBoundingBox. */
    public StsBoundingBox unrotatedBoundingBox;
    /** unrotated box whose range encompasses the unrotatedBoundingBox. Coordinates are "rounded" to even values. This is the box which is displayed. */
    public StsDisplayBoundingBox displayBoundingBox;
    /** defines grid for model grids/zones */
    public StsGridDefinition gridDefinition = null;
    /** Default Units */
    private byte timeUnits = StsParameters.TIME_MSECOND;
    private byte depthUnits = StsParameters.DIST_FEET;
    private byte xyUnits = StsParameters.DIST_FEET;
    /** Current model time - set by user. */
    private long projectTime;
    /** Model time range - set by user. */
    private long projectTimeDuration = 0;
    /** Default Colors */
    private StsColor defaultWellColor = new StsColor(StsColor.RED);
    private StsColor defaultSensorColor = new StsColor(StsColor.RED);
    /** Minimum depth; substituted into boundingBoxes when switched to depth */
    private float depthMin = StsParameters.largeFloat;
    /** Maximum depth; substituted into boundingBoxes when switched to depth */
    private float depthMax = -StsParameters.largeFloat;
    /** Depth increment; substituted into boundingBoxes when switched to depth */
    private float depthInc = 0.0f;
    /** Minimum time; substituted into boundingBoxes when switched to depth */
    private float timeMin = StsParameters.largeFloat;
    /** Maximum time; substituted into boundingBoxes when switched to depth */
    private float timeMax = -StsParameters.largeFloat;
    /** Time increment; substituted into boundingBoxes when switched to depth */
    private float timeInc;
    /** Depth increment locked: seismic loaded in depth; don't change when additional data loaded */
    public boolean depthIncLocked = false;
    /** Capture images in high (bmp) or low (jpg) resolution */
    public boolean highResolution = false;
    /** Time increment locked: seismic loaded in time; don't change when additional data loaded */
    public boolean timeIncLocked = false;
    /** defines cropped volume for project */
    public StsCropVolume cropVolume;
    /** domains which project supports: time, approx-depth, depth, time & depth, or approx-depth & depth */
    public byte zDomainSupported = TD_NONE;
    /** domain in which model was built if and when built */
    public byte modelZDomain = TD_NONE;
    /** current project domain (time or depth). */
    public byte zDomain = TD_DEPTH;
    /** Distance limit between project bounds and new data being added. If exceeded, probably data error */
    private double distanceLimit = 100000f;

    public float gridFrequency = AUTO;
    public float zGridFrequency = AUTO;
    public float labelFrequency = AUTO;

    public StsSeismicVelocityModel velocityModel = null;

    private float mapGenericNull = 1.e30f; // default to zmap null value
    private float logNull = -999;  // This is the value we expect in the file, not our internal null representation
    private float minZInc = 1.0f;

    private boolean displayContours = false;
    private float contourInterval = 10.0f;

    private boolean isVisible = true;
    private boolean showIntersection = true;
    private boolean showGrid = true;
    private boolean show3dGrid = false;
    private boolean show2dGrid = false;
    private boolean cursorPositionByMouse = true;

    // Time series plot defaults
    private boolean showTimeSeriesGrid = true;
    // Number of frames to leave accent up
    private StsColor timeSeriesBkgdColor = new StsColor(StsColor.BLACK);
    private StsColor timeSeriesFrgdColor = new StsColor(StsColor.LIGHT_GREY);
    private StsColor timeSeriesGridColor = new StsColor(StsColor.DARK_GREY);

    /** Show the labels in graphics environment */
    private boolean showLabels = false;
    private StsColor backgroundColor = new StsColor(StsColor.BLACK);
    private StsColor foregroundColor = new StsColor(StsColor.WHITE);
    private StsColor gridColor = new StsColor(StsColor.GRAY);
    private StsColor timingColor = new StsColor(StsColor.GRAY);
    private StsColor cogColor = new StsColor(StsColor.CYAN);
    private String labelFormatString = defaultLabelFormatString;
    private String timeDateFormatString = defaultDateFormatString + " " + defaultTimeFormatString;
    private String timeFormatString = defaultTimeFormatString;
    private String dateFormatString = defaultDateFormatString;
    private int dateOrder = CalendarParser.DD_MM_YY;

    private int timeUpdateRate = 1000; // Time action update rate when running time.
    private int realtimeBarAt = 500; // If more than 500 events in update, show progress bar

    // Project datums
    // private float timeDatum = 0.0f;
    // private float depthDatum = 0.0f;

    /** all 3d views use perspective; othewise orthogonal */
    private boolean isPerspective = true;
    public byte NOLIMIT = -1;
    private float zoomLimit = NOLIMIT; // Limit the zoom minimum
    /** Display Compass in 3D view*/
    private boolean displayCompass = true;

    transient StsColorComboBoxFieldBean backgroundColorBean;
    transient StsColorComboBoxFieldBean gridColorBean;
    transient StsColorComboBoxFieldBean timingColorBean;
    //  transient StsColorListFieldBean wellPanelColorBean;
    transient StsColorComboBoxFieldBean cogColorBean;
    transient StsBooleanFieldBean isPerspectiveBean;

    /** Manages memory for seismic classes */
    transient protected StsBlocksMemoryManager blocksMemoryManager = null;

    transient private DecimalFormat labelFormat = new DecimalFormat(defaultLabelFormatString);
    transient private SimpleDateFormat timeDateFormat = new SimpleDateFormat(defaultDateFormatString + " " + defaultTimeFormatString);
    transient private SimpleDateFormat dateFormat = new SimpleDateFormat(defaultDateFormatString);
    transient private SimpleDateFormat timeFormat = new SimpleDateFormat(defaultTimeFormatString);

    /**
     * these are stsClasses which are currently loaded (e.g. wells) before rotatedGrid angle is set (seismic or surfaces)
     * when and if the angle gets set, the coordinates of these unrotatedClass objects needs to be recomputed
     * in the rotated coordinate system.  Once this is done, these unrotatedClasses can be deleted from the list.
     */
    transient StsClass[] unrotatedClasses = null;

    transient private Properties userPreferences;

    transient private String rootDirname = null;
    transient private File rootDirectory = null;
    transient private File dataDirectory = new File(DATA_FOLDER);
    transient private File modelDirectory = new File(MODEL_FOLDER);
    transient private File binaryDirectory = new File(BINARY_FOLDER);
    transient private File archiveDirectory = new File(ARCHIVE_FOLDER);
    transient private File mediaDirectory = new File(MEDIA_FOLDER);
    transient private File exportModelsDirectory = new File(EXPORT_MODELS_FOLDER);
    transient private String name = "null";

    transient StsAsciiTokensFile directoryAliasesFile = null;

    transient StsComboBoxFieldBean defaultToolbarTimeDepthBean = null;

    transient StsTextureSurfaceFace basemapTextureObject = null;

	private transient float rowGridSize;
	private transient float colGridSize;
	private transient int nRows;
	private transient int nCols;
	transient float cellSize = 1.0f;

    static private String[] projectDirectories = null;
    static private String[] projectDatabases = null;
    static private String EOL = StsParameters.EOL;

    transient private StsObjectPanel objectPanel = null;

    static final boolean debug = true;

    /** Number of XY Grid increments */
    static public int approxNumXYGridIncrements = 20;

    /** Number of Z Grid increments */
    static public int approxNumZGridIncrements = 200;

    /** max number of incs allowed on vertical scale (time or depth) */
    static public int nZIncsMax = 1000;
    // file fields
    /** Default binary files folder name */
    static public final String BINARY_FOLDER = "BinaryFiles";

    /** Default model file folder name */
    static public final String MODEL_FOLDER = "modelDBFiles";

    /** Default archive file folder name */
    static public final String ARCHIVE_FOLDER = "archiveFiles";

    /** Default media file folder name */
    static public final String MEDIA_FOLDER = "mediaFiles";

    /** Default data files folder name */
    static public final String DATA_FOLDER = ".";
    static public final String EXPORT_MODELS_FOLDER = "exportedModelFiles";

    /** maximum allowable difference between rotation angles of objects in project */
    static public float angleTol = 0.5f;

    /** List of possible background colors */
    static public final StsColor[] backgroundColors = new StsColor[]
            {
                    StsColor.BLACK, StsColor.WHITE, StsColor.GRAY, StsColor.LIGHT_GRAY, StsColor.DARK_GRAY, StsColor.RED, StsColor.BLUE, StsColor.GREEN
            };
    /** foreground colors typically for text and ticks which contrasts with background color */
    static public final StsColor[] foregroundColors = new StsColor[]
            {
                    StsColor.WHITE, StsColor.BLACK, StsColor.WHITE, StsColor.BLACK, StsColor.WHITE, StsColor.WHITE, StsColor.WHITE, StsColor.WHITE
            };
    /** List of possible center-of-gravity colors */
    static public final StsColor[] cogColors = new StsColor[]
            {
                    StsColor.CYAN, StsColor.WHITE, StsColor.GRAY, StsColor.LIGHT_GRAY, StsColor.DARK_GRAY, StsColor.MAGENTA, StsColor.BLACK
            };

    static public final String defaultLabelFormatString = "#,###.#";
    static public final String defaultTimeFormatString = "HH:mm:ss.SSS";
    static public final String defaultDateFormatString = "dd/MM/yy";

    /** List of display beans for StsProject */
    transient protected StsFieldBean[] displayFields;

    /** List of property fields for StsProject */
    transient protected StsFieldBean[] propertyFields;

    /** explicitly-named propertyField used frequently */
    transient protected StsStringFieldBean zDomainSupportedBean = new StsStringFieldBean(this, "zDomainSupported", false, "Supported Domains: ");

    /** List of default fields for StsProject.  Currently none. */
    transient protected StsFieldBean[] defaultFields;

    static public final byte TD_DEPTH = StsParameters.TD_DEPTH; // 0
    static public final byte TD_TIME = StsParameters.TD_TIME; // 1
    static public final byte TD_TIME_DEPTH = StsParameters.TD_TIME_DEPTH; // 2
    static public final byte TD_APPROX_DEPTH_AND_DEPTH = StsParameters.TD_APPROX_DEPTH_AND_DEPTH; // 4
    static public final byte TD_NONE = StsParameters.TD_NONE; // 5

    static final byte defaultDistUnits = StsParameters.DIST_FEET;
    static final byte defaultTimeUnits = StsParameters.TIME_MSECOND;

	static final int maxGridSize = 1000;

    static final long serialVersionUID = 5811416831478613787L;
    public static final String PreStackSeismicDirectory = "PreStackSeismicDirectory";

    private transient StsAsciiTokensFile directoryTableFile;

    /** Default Project constructor. Uses user.dir to construct the Project object */
    public StsProject()
    {
        initializeProjectDirectory(System.getProperty("user.dir") + File.separator + "S2SCache" + File.separator);
    }

    public StsProject(boolean persistent)
    {
        super(persistent);
        initializeProjectDirectory(System.getProperty("user.dir") + File.separator + "S2SCache" + File.separator);
    }

    /**
     * Project constructor. Defines a default database name, directory and bounding box and initializes the graphics colors
     *
     * @parameters projectDirname the directory name where the database will be stored
     */
    public StsProject(String projectDirname)
    {
        super(false);
        initializeProjectDirectory(projectDirname);
    }

    private void initializeProjectDirectory(String projectDirname)
    {
        constructDefaultBoundingBoxes();
        File rootDir = new File(projectDirname);
        if(!rootDir.exists())
        {
            rootDir.mkdir();
        }
        setRootDirectory(rootDir);
        //TODO this does not seem to work for the directory; lastMod time is not changed
        //TODO we probably need to save project list file with two tokens, directory name and time modified
        //TODO time modified would be set when a db is saved in the project
        rootDir.setLastModified(System.currentTimeMillis());
        initializeColors();
    }

    public boolean initialize(StsModel model)
    {
        rootDirname = model.getDatabase().getProjectDirectory();
        if(!Main.isJarDB)
            rootDirectory = new File(rootDirname);
        else
            rootDirectory = new File(System.getProperty("user.home") + File.separator + "S2SCache" + File.separator);
        //        setColors(model);
        //constructDefaultBoundingBoxes();
        setLabelFormatString(labelFormatString);
        //        setZDomain(currentModel.getByteProperty("zDomain"));
        setIsDepth(zDomain == TD_DEPTH);

        return resetBoundingBoxesZRanges();

        //        supportedZDomainBean.setValueObject(StsParameters.TD_ALL_STRINGS[supportedZDomain]);
    }

    private StsColor getColorProperty(StsModel model, String key)
    {
        String colorString = model.properties.getProperty(key);
        return StsColor.colorFromString(colorString);
    }

    public void setLabelFormatString(String string)
    {
        if(string.equals(labelFormatString)) return;
        labelFormatString = new String(string);
        dbFieldChanged("labelFormatString", labelFormatString);
        labelFormat = new DecimalFormat(labelFormatString);
    }

    public void setDateOrderString(String dateOrderString)
    {
        dateOrder = CalendarParser.getOrderFromString(dateOrderString);
        switch(dateOrder)
        {
            case CalendarParser.DD_MM_YY:
                setDateFormatString("dd/MM/yy");
                break;
            case CalendarParser.DD_YY_MM:
                setDateFormatString("dd/yy/MM");
                break;
            case CalendarParser.MM_DD_YY:
                setDateFormatString("MM/dd/yy");
                break;
            case CalendarParser.MM_YY_DD:
                setDateFormatString("MM/yy/dd");
                break;
            case CalendarParser.YY_DD_MM:
                setDateFormatString("yy/dd/MM");
                break;
            case CalendarParser.YY_MM_DD:
                setDateFormatString("yy/MM/dd");
                break;
        }
        doObjectPanelChanged();
    }

    public String getDateOrderString()
    {
        return CalendarParser.getStringFromDateOrder(dateOrder);
    }

    public void setTimeFormatString(String string)
    {
        try
        {
            if(string.equals(timeFormatString)) return;
            timeFormat = new SimpleDateFormat(string);
            timeFormatString = new String(string);
            dbFieldChanged("timeFormatString", timeFormatString);
            resetTimeDateFormat();
        }
        catch(Exception e)
        {
            StsException.outputException("Error converting input string to time format", e, StsException.WARNING);
            timeFormat = new SimpleDateFormat(timeFormatString);
        }
    }

    public void setDateFormatString(String string)
    {
        try
        {
            if(string.equals(dateFormatString)) return;
            dateFormat = new SimpleDateFormat(string);
            dateFormatString = new String(string);
            dbFieldChanged("dateFormatString", dateFormatString);
            resetTimeDateFormat();
        }
        catch(Exception e)
        {
            StsException.outputException("Error converting input string to date format", e, StsException.WARNING);
            dateFormat = new SimpleDateFormat(dateFormatString);
        }
    }

    public void resetTimeDateFormat()
    {
        if(currentModel == null || currentModel.win3d == null) return;
        StsTimeActionToolbar timeActionToolbar = currentModel.win3d.getTimeActionToolbar();
        if(timeActionToolbar == null) return;
        timeActionToolbar.resetFormat();
        currentModel.resetTimeViews();
        currentModel.viewObjectRepaint(this, this);
    }

    public String getTimeDateFormatString()
    {
        return dateFormatString + " " + timeFormatString;
    }

    public String getDateFromLong(long time)
    {
        return getDateFromString(getTimeString(time));
    }

    public static String getDateFromString(String dateTime)
    {
        int firstColon = dateTime.indexOf(":");
        int lastSpace = dateTime.lastIndexOf(" ");
        if(firstColon > lastSpace)
            return dateTime.substring(0, lastSpace);
        else
        {
            String subString = dateTime.substring(0, lastSpace);
            lastSpace = subString.lastIndexOf(" ");
            return subString.substring(0, lastSpace);
        }
    }

    public int getDateOrder()
    {
        return dateOrder;
    }

    public boolean validateTimeDateString(String test)
    {
        try
        {
            getTimeDateFormat().parse(test);
            return true;
        }
        catch(Exception ex)
        {
            return false;
        }
    }

    public SimpleDateFormat getTimeDateFormat()
    {
        timeDateFormat = new SimpleDateFormat(dateFormatString + " " + timeFormatString);
        return timeDateFormat;
    }

    public SimpleDateFormat getTimeDateFormatForOutput()
    {
        timeDateFormat = new SimpleDateFormat(dateFormatString + ", " + timeFormatString);
        return timeDateFormat;
    }

    public SimpleDateFormat getDateFormat()
    {
        dateFormat = new SimpleDateFormat(dateFormatString);
        return dateFormat;
    }

    public SimpleDateFormat getTimeFormat()
    {
        timeFormat = new SimpleDateFormat(timeFormatString);
        return timeFormat;
    }

    public float getZoomLimit()
    {
        return zoomLimit;
    }

    public void setZoomLimit(int zLimit)
    {
        zoomLimit = zLimit;
        dbFieldChanged("zoomLimit", zoomLimit);
    }

    /** Set the Project Grid Frequency String */
    public void setZoomLimitString(String limit)
    {
        //		setDisplayField("zoomLimitAsString", limit);
        if(limit.equalsIgnoreCase("No") || limit.equalsIgnoreCase("No Limit") || limit.equalsIgnoreCase("None"))
        {
            zoomLimit = NOLIMIT;
        }
        else
        {
            zoomLimit = (new Float(limit)).floatValue();
        }
        //		dbFieldChanged("zoomLimit",zoomLimit);
        return;
    }

    /**
     * Get Project Label Frequency String
     *
     * @return label frequency in model units
     */
    public String getZoomLimitString()
    {
        if(zoomLimit == NOLIMIT)
        {
            return new String("No Limit");
        }
        else
        {
            return new Float(zoomLimit).toString();
        }
    }

    public double getDistanceLimit()
    {
        return distanceLimit;
    }

    public void setDistanceLimit(double limit)
    {
        if(limit == distanceLimit) return;
        distanceLimit = limit;
        dbFieldChanged("distanceLimit", distanceLimit);
    }

    public int getRealtimeBarAt()
    {
        return realtimeBarAt;
    }

    public void setRealtimeBarAt(int numEvents)
    {
        if(numEvents == realtimeBarAt) return;
        realtimeBarAt = numEvents;
        dbFieldChanged("realtimeBarAt", realtimeBarAt);
    }

    public int getTimeUpdateRate()
    {
        return timeUpdateRate;
    }

    public void setTimeUpdateRate(int rate)
    {
        if(rate == timeUpdateRate) return;
        timeUpdateRate = rate;
        dbFieldChanged("timeUpdateRate", timeUpdateRate);
    }

    public DecimalFormat getLabelFormat()
    {
        if(labelFormat == null)
            labelFormat = new DecimalFormat(labelFormatString);
        return labelFormat;
    }

    /** Initialize the background, grid and timing line colors to black, light gray, and gray respectively */
    private void initializeColors()
    {
        backgroundColor = new StsColor(StsColor.BLACK);
        gridColor = new StsColor(StsColor.LIGHT_GRAY);
        timingColor = new StsColor(StsColor.GRAY);
        cogColor = new StsColor(StsColor.CYAN);
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;

        backgroundColorBean = new StsColorComboBoxFieldBean(this, "backgroundColor", "Background Color: ", backgroundColors);
        gridColorBean = new StsColorComboBoxFieldBean(this, "gridColor", "Grid Color: ", backgroundColors);
        timingColorBean = new StsColorComboBoxFieldBean(this, "timingColor", "Timing Color: ", backgroundColors);
        //wellPanelColorBean = new StsColorListFieldBean(this, "wellPanelColor", "Well panel Color:", backgroundColors);
        cogColorBean = new StsColorComboBoxFieldBean(this, "cogColor", "COG Color: ", cogColors);
        isPerspectiveBean = new StsBooleanFieldBean(this, "isPerspective", "Perspective View");
        displayFields = new StsFieldBean[]
                {
                        new StsBooleanFieldBean(this, "isVisible", "Show Bounding Box"),
                        new StsBooleanFieldBean(this, "displayCompass", "Display Compass"),
                        isPerspectiveBean,
                        new StsBooleanFieldBean(this, "showGrid", "Show Floor Grid"),
                        new StsBooleanFieldBean(this, "showIntersection", "Cursor Intersection"),
                        new StsBooleanFieldBean(this, "cursorPositionByMouse", "Position Cursor with Mouse"),                        
                        backgroundColorBean,
                        gridColorBean,
                        timingColorBean,
                        new StsBooleanFieldBean(this, "enableTsGrid", "Show Grid on Series Plot"),
                        new StsColorComboBoxFieldBean(this, "timeSeriesBkgdColor", "Series Plot Background Color:", StsColor.bkgdColors),
                        new StsColorComboBoxFieldBean(this, "timeSeriesFrgdColor", "Series Plot Foreground Color:", StsColor.bkgdColors),
                        new StsColorComboBoxFieldBean(this, "timeSeriesGridColor", "Series Plot Grid Color:", StsColor.bkgdColors),
                        //Bean,
                        new StsBooleanFieldBean(this, "showLabels", "Show Labels on Cursor"),
                        new StsBooleanFieldBean(this, "show3dGrid", "Show Grid in 3D Views"),
                        new StsBooleanFieldBean(this, "show2dGrid", "Show Grid in 2D Views"),
                        cogColorBean,
                        new StsStringFieldBean(this, "labelFormatString", "Label Format: "),
                        new StsStringFieldBean(this, "labelFrequencyString", "Label Frequency: "),
                        new StsStringFieldBean(this, "gridFrequencyString", "XY Grid Frequency: "),
                        new StsStringFieldBean(this, "zGridFrequencyString", "Z Grid Frequency: "),
                        new StsComboBoxFieldBean(this, "dateOrderString", "Date Order: ", CalendarParser.dateOrderStrings),
                        new StsStringFieldBean(this, "timeFormatString", "Time Format: "),
                        new StsStringFieldBean(this, "dateFormatString", "Date Format: "),
                        new StsIntFieldBean(this, "timeUpdateRate", 100, 600000, "Time Update Rate(ms): ", true),
                        new StsIntFieldBean(this, "realtimeBarAt", 100, 100000, "Progress Bar At (nEvents): ", true),
                        new StsDoubleFieldBean(this, "distanceLimit", 1000, 1000000, "Distance Limit (ft/m): ", true),
                        new StsStringFieldBean(this, "zoomLimitString", "Zoom Limit(ft/m): "),
                        new StsBooleanFieldBean(this, "highResolution", "Capture High Res Images"),
                        new StsFloatFieldBean(this, "mapGenericNull", true, "Surface Null Z: "),
                        new StsFloatFieldBean(this, "logNull", "Log Null Value: "),
                        new StsFloatFieldBean(this, "minZInc", "Min time or depth inc: "),
                        new StsBooleanFieldBean(this, "displayContours", "Display Contours"),
                        new StsFloatFieldBean(this, "contourInterval", "Contour Interval: ")
                };

        return displayFields;
    }

      public StsFieldBean[] getPropertyFields()
      {
          if(propertyFields == null)
          {
			  StsFloatFieldBean cellSizeBean = new StsFloatFieldBean(this, "cellSize", 0.0f, 1000.0f, "Cell Size:", true);
    		  cellSizeBean.fixStep(1.0f);
    		  cellSizeBean.setToolTipText("Specify the cell size within the stimulated volume.");
              propertyFields = new StsFieldBean[]
                      {
                              //new StsFloatFieldBean(this, "timeDatum", true, "Time Datum: "),
                              //new StsFloatFieldBean(this, "depthDatum", true, "Depth Datum: "),
                              new StsFloatFieldBean(this, "displayXMin", false, "X Minimum: "),
                              new StsFloatFieldBean(this, "displayXMax", false, "X Maximum: "),
                              new StsFloatFieldBean(this, "displayYMin", false, "Y Minimum: "),
                              new StsFloatFieldBean(this, "displayYMax", false, "Y Maximum: "),
                              new StsFloatFieldBean(this, "angle", false, "XY Rot Angle: "),
                              //		new StsFloatFieldBean(this, "zStart", "Start Z", false),
                              new StsDoubleFieldBean(this, "xOrigin", false, "X Origin: "),
                              new StsDoubleFieldBean(this, "yOrigin", false, "Y Origin: "),
                              new StsFloatFieldBean(this, "userDepthMin", true, "Depth Min: "),
                              new StsFloatFieldBean(this, "userDepthMax", true, "Depth Max: "),
                              new StsFloatFieldBean(this, "depthInc", false, "Depth Inc: "),
                              //new StsButtonFieldBean("Reset","Reset the bounding box to user supplied limits.",this,"resetProjectBounds"),
                              new StsFloatFieldBean(this, "timeMin", false, "Time Min: "),
                              new StsFloatFieldBean(this, "timeMax", false, "Time Max: "),
                              new StsFloatFieldBean(this, "timeInc", false, "Time Inc: "),

                              new StsStringFieldBean(this, "timeUnitString", false, "Time Units: "),
                              new StsStringFieldBean(this, "depthUnitString", false, "Depth Units: "),
                              new StsStringFieldBean(this, "xyUnitString", false, "Distance Units: "),
                              //		zDomainBean,
                  			  zDomainSupportedBean,
							  cellSizeBean,
                              new StsFloatFieldBean(this, "gridDX", false, "Grid dX: "),
                              new StsFloatFieldBean(this, "gridDY", false, "Grid dY: "),
                              new StsFloatFieldBean(this, "gridZ", true, "Grid Z: ", true)
                              //        new StsFloatFieldBean(this, "cultureZ", true, "Culture Z", true)
                      };
          }
          return propertyFields;
      }

      public StsFieldBean[] getDefaultFields() { return defaultFields; }

      /** Set the current Model related to this Project
       * @parameters model the current model
       */
  //    public void setModel(StsModel model) { this.model = model; }

      /** Create the default bounding box */
      private void constructDefaultBoundingBoxes()
      {
          rotatedBoundingBox = new StsRotatedGridBoundingBox(false, "Project.rotatedGridBoundingBox");
          rotatedBoundingBox.initialize(-50f, 50f, 1.0f, -50f, 50f, 1.0f, 0f, 100f, 1.0f, false);

          unrotatedBoundingBox = new StsBoundingBox(false, "Project.unrotatedBoundingBox");
          unrotatedBoundingBox.initialize(-50f, 50f, -50f, 50f, 0f, 100f, false);
          displayBoundingBox = new StsDisplayBoundingBox(false, "Project.displayBoundingBox");
          displayBoundingBox.initialize(-50f, 50f, -50f, 50f, 0f, 100f, false);
          depthMin = 0f;
          depthMax = 100f;
          depthInc = 1.0f;
          setZDomain(TD_DEPTH);
          initializeBoundingBoxes();
          cropVolume = new StsCropVolume(false);
      }

      /**
       * after db is defined, we need to set non-persistent boundingBoxes to persistent
       * so they will be written to db.
       */
      public void dbInitialize(StsModel model)
      {
          addToModel();
          if(rotatedBoundingBox != null)
          {
              rotatedBoundingBox.addToModel();
              dbFieldChanged("rotatedBoundingBox", rotatedBoundingBox);
          }
          if(unrotatedBoundingBox != null)
          {
              unrotatedBoundingBox.addToModel();
              dbFieldChanged("unrotatedBoundingBox", unrotatedBoundingBox);
          }
          if(displayBoundingBox != null)
          {
              displayBoundingBox.addToModel();
              dbFieldChanged("displayBoundingBox", displayBoundingBox);
          }
          if(cropVolume != null)
          {
              cropVolume.addToModel();
              dbFieldChanged("cropVolume", cropVolume);
          }
 //         initializeZDomainAndRebuild();
      }

      /** call when any displayed project parameters have changed to update object panel */
      public void objectPanelChanged()
      {
          StsToolkit.runLaterOnEventThread(new Runnable()
          {
              public void run() { doObjectPanelChanged(); }
          });
      }

      private void doObjectPanelChanged()
      {
          objectPanel = getObjectPanel();
          objectPanel.refreshProperties();
      }

      /**
       * Get the project directories from the S2S user projects file. These are the directories to be displayed in the
       * File drop down.
       *
       * @returns projectDirectories a string array of project directories
       */
      public String[] getProjectDirectories()
      {
          if(projectDirectories == null)
          {
              String directory = getUserPropertiesDirectory();
              String filename = "s2s.user.projects";
              StsMessageFiles.infoMessage("Preferences will be stored in directory: " + directory);
              projectDirectories = readPathnamesFromFile(directory, filename);
          }
          return projectDirectories;
      }

      static public String getUserPropertiesDirectory()
      {
          return System.getProperty("user.home") + File.separator + "S2SCache" + File.separator;
      }

      /**
       * Add a project directory to the S2S user projects file. These directories will be displayed in the
       * File drop down.
       *
       * @parameters projectName the directory name
       */
      static public void addProjectDirectory(String projectName)
      {
          int nProjects = projectDirectories.length;
          for(int n = 0; n < nProjects; n++)
          {
              if(projectName.equals(projectDirectories[n]))
              {
                  return;
              }
          }

          projectDirectories = (String[])StsMath.arrayInsertElementBefore(projectDirectories, projectName, 0);
          String directory = getUserPropertiesDirectory();
          String filename = "s2s.user.projects";
          StsFile.checkCreateDirectory(directory);
          StsFile file = StsFile.constructor(directory, filename);
          if(nProjects > 10)
              projectDirectories = (String[])StsMath.trimArray(projectDirectories, 10);
          file.writeStringsToFile(projectDirectories);
      }

      /**
       * Get all the current user databases
       *
       * @returns projectDatabases an array of valid databases read from the user database file
       */
      static public String[] getDatabases()
      {
          if(projectDatabases == null)
          {
              String directory = getUserPropertiesDirectory();
              String filename = "s2s.user.databases";
              projectDatabases = readPathnamesFromFile(directory, filename);
          }
          return projectDatabases;
      }

      /**
       * Initialize the database information and add the database to the list of user databases
       *
       * @parameters dbFile the new Sts database file
       */
      public void setDatabaseInfo(StsDBFile dbFile)
      {
          try
          {
              //String modelDbPathname = dbFile.getURLDirectory();
              //rootDirname = StsFile.getDirectoryFromPathname(modelDbPathname);
              //setRootDirectory(new File(rootDirname));
              String filename = dbFile.getFilename();
              name = filename.substring(filename.indexOf("db.")+3, filename.length());
              addDbPathname(dbFile.getURLPathname());
          }
          catch(Exception e)
          {
              StsException.outputException("StsProject.setDatabase() failed.",
                      e, StsException.WARNING);
          }
      }

      /**
       * Add a database pathname to the S2S user databases file.
       * The list is ordered from the most recent to the oldest.
       *
       * @parameters dbPathname pathname of the database
       */
      static public void addDbPathname(String dbPathname)
      {
          try
          {
              if(dbPathname == null) return;
              int nProjectDBs = 0;
              if(projectDatabases != null)
              {
                  nProjectDBs = projectDatabases.length;
                  for(int n = 0; n < nProjectDBs; n++)
                      if(dbPathname.equals(projectDatabases[n])) return;
              }

              projectDatabases = (String[])StsMath.arrayInsertElementBefore(projectDatabases, dbPathname, 0);
              String directory = getUserPropertiesDirectory();
              String filename = "s2s.user.databases";
              if(!StsFile.checkCreateDirectory(directory)) return;
              StsFile file = StsFile.constructor(directory, filename);
              if(projectDatabases.length > 20)
                  projectDatabases = (String[])StsMath.trimArray(projectDatabases, 20);
              file.writeStringsToFile(projectDatabases);
          }
          catch(Exception e)
          {
              StsException.outputException("StsProject.addDbPathname() failed.",
                      e, StsException.WARNING);
          }
      }

      /**
       * Reads pathnames from the specified file. Validates that the line read is a file.
       *
       * @parameters directory directory name
       * @parameters filename file name
       * @returns pathnames an array of valid paths read from the file
       */
      static public String[] readPathnamesFromFile(String directory, String filename)
      {
          String[] pathnames = new String[0];
          String pathname;
          boolean rewrite = false;
          StsFile.checkCreateDirectory(directory);

          StsFile file = StsFile.constructor(directory, filename);
          StsAsciiFile asciiFile = new StsAsciiFile(file);
          if(!asciiFile.openRead()) return pathnames;

          try
          {
              File f;

              while((pathname = asciiFile.readLine()) != null)
              {
                  if(pathname.startsWith("file:"))
                  {
                      URL url = new URL(pathname);
                      f = new File(url.getFile());
                  }
                  else
                  {
                      f = new File(pathname);

                  }
                  if(f == null || !f.exists())
                  {
                      rewrite = true;
                      continue;
                  }
                  pathnames = (String[])StsMath.arrayAddElement(pathnames, pathname);
              }
              if(rewrite)
              {
                  file.writeStringsToFile(pathnames);
              }
              asciiFile.close();
              return pathnames;
          }
          catch(Exception e)
          {
              StsException.systemError("Failed to read " + directory + filename);
              if(asciiFile != null)
              {
                  asciiFile.close();
              }
              return null;
          }
      }

      /*
        OutputStream os = new OutputStream();
=======
    /*
        OutputStream os = new OutputStream();
>>>>>>> 1.314
        file.
        PersistentStrings persistentStrings = PersistentStrings.read(filename);
        ArrayList projectDirectories = persistentStrings.getStringsArrayList();
       rootDirname = StsSelectProjectDirectoryDialog.getProjectDirectory("Select project directory", projectDirectories);
        rootDirectory = new File(rootDirname);
        projectDirectories.remove(rootDirname);
        projectDirectories.add(0, rootDirname);
        persistentStrings.write(filename);
       }
       catch(Exception e)
       {
        StsException.outputException("StsProject.getUserPreferences() failed.",
       e, StsException.WARNING);
       }
       }
       */
    /*
       public void loadUserPreferences()
       {
        try
        {
         String filename = System.getProperty("java.home") + File.separator + "s2s.user.projects";
         PersistentStrings persistentStrings = PersistentStrings.read(filename);
         ArrayList projectDirectories = persistentStrings.getStringsArrayList();
       rootDirname = StsSelectProjectDirectoryDialog.getProjectDirectory("Select project directory", projectDirectories);
         rootDirectory = new File(rootDirname);
         projectDirectories.remove(rootDirname);
         projectDirectories.add(0, rootDirname);
         persistentStrings.write(filename);
        }
        catch(Exception e)
        {
         StsException.outputException("StsProject.getUserPreferences() failed.",
        e, StsException.WARNING);
        }
       }
       */

    /**
     * Add a user preference to the database
     *
     * @parameters key the name to refer to the user preference
     * @parameters value the current value for the key
     */
    public void addUserPreference(String key, String value)
    {
        userPreferences.setProperty(key, value);
    }

    /**
     * Get a user preference from the database
     *
     * @parameters key the name to refer to the user preference
     * @parameters value the current value for the key
     */
    public String getUserPreference(String key)
    {
        return userPreferences.getProperty(key);
    }

    /**
     * Has the bounding box origin been explicitly set
     *
     * @return true if set
     */
    public boolean isOriginSet()
    {
        return rotatedBoundingBox.originSet;
    }

    public boolean checkSetOrigin(double xOrigin, double yOrigin)
    {
        if(unrotatedBoundingBox == null)
        {
            return false;
        }
        if(unrotatedBoundingBox.originSet)
        {
            return false; // origin already set: return false
        }
        initializeOriginAndRange(xOrigin, yOrigin);
        return true;
    }

    /** set origin of bounding boxes and reinitialize the boundingBox limits */
    private void initializeOriginAndRange(double xOrigin, double yOrigin)
    {
        unrotatedBoundingBox.initializeOriginAndRange(xOrigin, yOrigin);
        rotatedBoundingBox.initializeOriginAndRange(xOrigin, yOrigin);
        depthMin = StsParameters.largeFloat;
        depthMax = -StsParameters.largeFloat;
        timeMin = StsParameters.largeFloat;
        timeMax = -StsParameters.largeFloat;
    }

    /**
     * if angle is already set, and new angle is not the same, return false indicated this object can't be loaded;
     * if angle is not already set and new angle is
     */
    public boolean checkSetAngle(float angle)
    {

        if(rotatedBoundingBox == null) return false;

        boolean angleSet = rotatedBoundingBox.getAngleSet();
        float currentAngle = rotatedBoundingBox.angle;
        if(angleSet)
        {
            if(!StsMath.sameAsTol(currentAngle, angle, angleTol))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING,
                        "Grid angle: " + angle + " differs from current project angle: " + currentAngle + "\n" +
                                "by more than " + angleTol + " degrees.\n" +
                                "Can't load this object.");
                return false;
            }
        }
        else
        {
            rotatedBoundingBox.setAngle(angle);

        }
        if(unrotatedClasses != null && currentAngle != angle)
        {
            for(int n = 0; n < unrotatedClasses.length; n++)
            {
                unrotatedClasses[n].projectRotationAngleChanged();
            }
        }
        unrotatedClasses = null;
        return true;
    }

    /**
     * Initialize the rotated bounding box by setting the origin and rotation angle
     *
     * @parameters xOrigin actual x origin
     * @parameters yOrigin actual y origin
     * @parameters angle the angle off true North
     */
    public void setOriginAndAngle(double xOrigin, double yOrigin, float angle)
    {
        rotatedBoundingBox.checkSetOriginAndAngle(xOrigin, yOrigin, angle);
    }

    public void setOrigin(double xOrigin, double yOrigin)
    {
        rotatedBoundingBox.setOrigin(xOrigin, yOrigin);
    }
    /**
     * Compute the relative origin of the rotated bounding box
     *
     * @return xOrigin, yOrigin
     */
    public float[] computeRelativeOrigin(double xOrigin, double yOrigin)
    {
        if(rotatedBoundingBox == null)
        {
            return null;
        }
        else
        {
            return rotatedBoundingBox.computeRelativeOrigin(xOrigin, yOrigin);
        }
    }

    public boolean addToProject(StsCultureObjectSet2D cSet)
    {
        try
        {
            byte zDomainSupported = cSet.getZDomainSupported();
            setZDomainSupported(TD_TIME_DEPTH); // Always both for culture since it defaults to z of current slice.
            float newTimeMin = (float) cSet.getMinTime();
            float newTimeMax = (float) cSet.getMaxTime();
            float newDepthMin = (float) cSet.getMinDepth();
            float newDepthMax = (float) cSet.getMaxDepth();

            double xOrigin = cSet.getXMax();
            double yOrigin = cSet.getYMax();

            checkSetOrigin(xOrigin, yOrigin);

            //boolean supportsTime = supportsTime(supportedZDomain);
            //boolean supportsDepth = supportsDepth(supportedZDomain);

            int nObjectsInSet = cSet.getCultureObjects().length;
            for(int i = 0; i < nObjectsInSet; i++)
            {
                StsCultureObject2D object = cSet.getCultureObjects()[i];
                for(int j = 0; j < object.getNumPoints(); j++)
                    unrotatedBoundingBox.addXY(object.getPointAt(j), this.getXOrigin(), this.getYOrigin());
            }

            setTimeMin(newTimeMin);
            setTimeMax(newTimeMax);

            double[] scale = StsMath.niceScale(newDepthMin, newDepthMax, 100, true);
            setDepthMin((float) scale[0]);
            setDepthMax((float) scale[1]);
            checkSetDepthInc((float) scale[2]);

            if(zDomain == TD_DEPTH)
            {
                unrotatedBoundingBox.zMin = depthMin;
                unrotatedBoundingBox.zMax = depthMax;
            }
            else if(zDomain == TD_TIME)
            {
                unrotatedBoundingBox.zMin = timeMin;
                unrotatedBoundingBox.zMax = timeMax;
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.addToProject(StsLine) failed.", e, StsException.WARNING);
            return false;
        }
    }

    public double[] niceScale(float zMin, float zMax)
    {
        return StsMath.niceScale(zMin, zMax, approxNumZGridIncrements, true);
    }

    public boolean addToProject(StsLine line)
    {
        try
        {
            byte zDomainSupported = line.getZDomainSupported();
            if(!checkSetZDomain(zDomainSupported)) return false;
            double xOrigin = line.getXOrigin();
            double yOrigin = line.getYOrigin();
            checkSetOrigin(xOrigin, yOrigin);
            StsObjectRefList lineVertices = line.getLineVertices();
            int nVertices = lineVertices.getSize();
            float newTimeMin = StsParameters.largeFloat;
            float newTimeMax = -StsParameters.largeFloat;
            float newDepthMin = StsParameters.largeFloat;
            float newDepthMax = -StsParameters.largeFloat;
            boolean supportsTime = supportsTime(zDomainSupported);
            boolean supportsDepth = supportsDepth(zDomainSupported);
            boolean loadAnyway = false;
            for(int n = 0; n < nVertices; n++)
            {
                StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
                StsPoint point = vertex.getPoint();
                if(!unrotatedBoundingBox.sanityCheck(point, xOrigin, yOrigin, distanceLimit) && !loadAnyway)
                {
                    if(!StsYesNoDialog.questionValue(currentModel.win3d, "Line (" + line.getName() + " has data out of range limit (" + distanceLimit + " " + this.getXyUnitString() + ").\n\n" +
                            "   Do you want to continue?\n"))
                        return false;
                    else
                        loadAnyway = true;
                }
                unrotatedBoundingBox.addXY(point, xOrigin, yOrigin);
                if(supportsTime(zDomainSupported))
                {
                    float time = point.getT();
                    newTimeMin = Math.min(newTimeMin, time);
                    newTimeMax = Math.max(newTimeMax, time);
                }
                if(supportsDepth(zDomainSupported))
                {
                    float depth = point.getZ();
                    newDepthMin = Math.min(newDepthMin, depth);
                    newDepthMax = Math.max(newDepthMax, depth);
                    //                    checkDepthRange();
                }
            }

            if(supportsTime)
            {
                setTimeMin(newTimeMin);
                setTimeMax(newTimeMax);
            }
            if(supportsDepth)
            {
                double[] scale = niceScale(newDepthMin, newDepthMax);

                setDepthMin((float) scale[0]);
                setDepthMax((float) scale[1]);
                checkSetDepthInc((float) scale[2]);
            }
            if(zDomain == TD_DEPTH)
            {
                unrotatedBoundingBox.zMin = depthMin;
                unrotatedBoundingBox.zMax = depthMax;
            }
            else if(zDomain == TD_TIME)
            {
                unrotatedBoundingBox.zMin = timeMin;
                unrotatedBoundingBox.zMax = timeMax;
            }
            //            objectPanelChanged();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.addToProject(StsLine) failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean addToProject(StsBoundingBox boundingBox, byte zdomain)
    {
        if(!checkSetZDomain(zdomain))
        {
            return false;
        }
        if(!unrotatedBoundingBox.sanityCheck(boundingBox, distanceLimit))
        {
            if(!StsYesNoDialog.questionValue(currentModel.win3d, boundingBox.getName() + " has data out of range limit (" + distanceLimit + " " + this.getXyUnitString() + ").\n\n" +
                    "   Do you want to continue?\n"))
                return false;
        }
        unrotatedBoundingBox.addUnrotatedBoundingBox(boundingBox);
        double[] scale = niceScale(boundingBox.zMin, boundingBox.zMax);
        if(zDomain == TD_TIME)
        {
            setTimeMin((float) scale[0]);
            setTimeMax((float) scale[1]);
        }
        else
        {
            setDepthMin((float) scale[0]);
            setDepthMax((float) scale[1]);
            checkSetDepthInc((float) scale[2]);
        }
        if(this.zDomain == TD_DEPTH)
        {
            unrotatedBoundingBox.zMin = depthMin;
            unrotatedBoundingBox.zMax = depthMax;
        }
        else if(this.zDomain == TD_TIME)
        {
            unrotatedBoundingBox.zMin = timeMin;
            unrotatedBoundingBox.zMax = timeMax;
        }
        return true;
    }

    public boolean addToProject(StsSensor sensor, byte zdomain)
    {
        if(!checkSetZDomain(zdomain))
        {
            return false;
        }
        //checkSetOrigin(sensor.getXOrigin(),sensor.getYOrigin());
        StsObjectRefList timeCurves = sensor.getTimeCurves();
        int nCurves = timeCurves.getSize();

        float newDepthMin = StsParameters.largeFloat;
        float newDepthMax = -StsParameters.largeFloat;
        try
        {
            if(sensor instanceof StsStaticSensor)
            {
                double xOrig = 0.0f;
                double yOrig = 0.0f;

                float depth = (float) ((StsStaticSensor) sensor).getZLoc();
                newDepthMin = depth;
                newDepthMax = depth;
                if(newDepthMin < 0)
                    newDepthMin = depth + 1.0f;
                else
                    newDepthMin = depth - 1.0f;

                if(newDepthMax < 0)
                    newDepthMax = depth - 1.0f;
                else
                    newDepthMax = depth + 1.0f;

                checkSetOrigin(xOrig, yOrig);
                StsPoint point = new StsPoint(((StsStaticSensor) sensor).getXLoc(), ((StsStaticSensor) sensor).getYLoc(), ((StsStaticSensor) sensor).getZLoc());
                if(!unrotatedBoundingBox.sanityCheck(((StsStaticSensor) sensor).getXLoc(), ((StsStaticSensor) sensor).getYLoc(), xOrig, yOrig, distanceLimit))
                {
                    if(!StsYesNoDialog.questionValue(currentModel.win3d, sensor.getName() + " has data out of range limit (" + distanceLimit + " " + this.getXyUnitString() + ").\n\n" +
                            "   Do you want to continue?\n"))
                        return false;
                }
                unrotatedBoundingBox.addXY(point, xOrig, yOrig);
                //rotatedBoundingBox.addPoint(new float[] {(float)sensor.getXLoc(), (float)sensor.getYLoc(), (float)sensor.getZLoc()});
            }
            else
            {
                float[] depthFloats = null;
                float[] xFloats = null;
                float[] yFloats = null;
                double xOrig = 0.0f;
                double yOrig = 0.0f;

                for(int i = 0; i < nCurves; i++)
                {
                    StsTimeCurve timeCurve = (StsTimeCurve) timeCurves.getElement(i);
                    if(timeCurve.getValueVector().getType() == StsLogVector.DEPTH)
                        depthFloats = timeCurve.getValuesVectorFloats();
                    else if(timeCurve.getValueVector().getType() == StsLogVector.X)
                    {
                        xFloats = timeCurve.getValuesVectorFloats();
                        xOrig = timeCurve.getValueVector().getOrigin();
                    }
                    else if(timeCurve.getValueVector().getType() == StsLogVector.Y)
                    {
                        yFloats = timeCurve.getValuesVectorFloats();
                        yOrig = timeCurve.getValueVector().getOrigin();
                    }
                }
                checkSetOrigin(xOrig, yOrig);
                if(xFloats != null) // Realtime will be initialized with no values.
                {
                    boolean loadAnyway = false;
                    for(int n = 0; n < xFloats.length; n++)
                    {
                        if(!unrotatedBoundingBox.sanityCheck(xFloats[n], yFloats[n], xOrig, yOrig, distanceLimit) && !loadAnyway)
                        {
                            if(!StsYesNoDialog.questionValue(currentModel.win3d, "Data Out of Range Limit of " + distanceLimit + " " + this.getXyUnitString() + ".\n\n" +
                                    "   Do you want to continue?\n"))
                                return false;
                            else
                                loadAnyway = true;
                        }
                        StsPoint point = new StsPoint(xFloats[n], yFloats[n], depthFloats[n]);
                        unrotatedBoundingBox.addXY(point, xOrig, yOrig);
                        newDepthMin = Math.min(newDepthMin, depthFloats[n]);
                        newDepthMax = Math.max(newDepthMax, depthFloats[n]);
                    }
                    if(newDepthMin == newDepthMax)
                    {
                        newDepthMax = newDepthMin + 1;
                        StsMessageFiles.infoMessage("Unable to fit 3D box around data, Z has no dimension. Adding 1 to the maximum");
                    }
                }
            }

            if(sensor instanceof StsDynamicSensor)
            {
                double[] scale = niceScale(newDepthMin, newDepthMax);
                /*
                          if(zdomain == TD_DEPTH)
                          {
                              if(scale[0] < depthMin) setDepthMin((float) scale[0]);
                              if(scale[1] > depthMax) setDepthMax((float) scale[1]);
                              if(scale[2] != depthInc) checkSetDepthInc((float) scale[2]);
                          }
                          else
                          {
                              if(scale[0] < timeMin) setTimeMin((float) scale[0]);
                              if(scale[0] < timeMin) setTimeMax((float) scale[1]);
                              if(scale[0] < timeMin) checkSetTimeInc((float) scale[2]);
                          }
              */

                if(zDomain == TD_TIME)
                {
                    setTimeMin((float) scale[0]);
                    setTimeMax((float) scale[1]);
                }
                else
                {
                    setDepthMin((float) scale[0]);
                    setDepthMax((float) scale[1]);
                    checkSetDepthInc((float) scale[2]);
                }
                if(this.zDomain == TD_DEPTH)
                {
                    unrotatedBoundingBox.zMin = depthMin;
                    unrotatedBoundingBox.zMax = depthMax;
                }
                else if(this.zDomain == TD_TIME)
                {
                    unrotatedBoundingBox.zMin = timeMin;
                    unrotatedBoundingBox.zMax = timeMax;
                }
            }
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.addToProject(StsPoint[]) failed.", e, StsException.WARNING);
            return false;
        }
    }

    static public boolean supportsTime(byte zDomain)
    {
        return zDomain == TD_TIME_DEPTH || zDomain == TD_TIME;
    }

    static public boolean supportsDepth(byte zDomain)
    {
        return zDomain == TD_TIME_DEPTH || zDomain == TD_DEPTH;
    }

    public boolean supportsTime()
    {
        return supportsTime(zDomainSupported);
    }

    public boolean supportsDepth()
    {
        return supportsDepth(zDomainSupported);
    }

    public boolean supportsZDomain(byte objectZDomain)
    {
        return zDomainSupported == TD_TIME_DEPTH || objectZDomain == zDomain;
    }

    public boolean canDisplayZDomain(byte objectZDomain)
    {
        return objectZDomain == TD_TIME_DEPTH || objectZDomain == zDomain;
    }

    public boolean hasVelocityModel()
    {
        if(velocityModel == null)
            return false;
        else
            return true;
    }

    /*
          public boolean addZToProject(StsSeismicBoundingBox box)
          {
              byte objectZDomain = StsParameters.getZDomainFromString(box.zDomain);
              if (objectZDomain == TD_DEPTH)
                  {
                      depthMin = rotatedBoundingBox.zMin;
                      depthMax = rotatedBoundingBox.zMax;
                      depthInc = rotatedBoundingBox.zInc;
      //				fieldChanged("depthIncLocked", true);
                  }
                  else if (zDomain == TD_TIME)
                  {
                      timeMin = rotatedBoundingBox.zMin;
                      timeMax = rotatedBoundingBox.zMax;
                      timeInc = rotatedBoundingBox.zInc;
      //				fieldChanged("timeIncLocked", true);
                  }
                  rangeChanged();
                  objectPanelChanged();
             return true;
          }
      */
    //TODO need to have a flag saying a time or depth congruent volume has been loaded and any changes to project z range need to be congruent with this
    public boolean addToProject(StsSeismicBoundingBox volume, boolean setAngle)
    {
        try
        {
            // check data

            if(volume.zMax <= volume.zMin)
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "Volume zMax " + volume.zMax + " <= zMin " + volume.zMin + "\nNot loading volume.");
                return false;
            }
            byte objectZDomain = StsParameters.getZDomainFromString(volume.zDomain);
            if(!checkSetZDomain(objectZDomain))
            {
                return false;
            }

            double volumeXOrigin = volume.xOrigin;
            double volumeYOrigin = volume.yOrigin;

            boolean originChanged = checkSetOrigin(volumeXOrigin, volumeYOrigin);
            if(setAngle)
            {
                boolean angleAlreadySet = rotatedBoundingBox.getAngleSet();
                if(!checkSetAngle(volume.angle)) return false;

                // If angle not already set, then the current rotatedBoundingBox is same as the
                // project.unrotatedBoundingBox. So if we have set the angle now, we start with
                // a new unrotatedBoundingbox and add this rotated object to it.
                if(!angleAlreadySet)
                {
                    rotatedBoundingBox.resetRange();
                }
            }
            float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY(volumeXOrigin, volumeYOrigin);
            volume.xMin += xy[0];
            volume.yMin += xy[1];
            volume.xMax += xy[0];
            volume.yMax += xy[1];
            volume.xOrigin = getXOrigin();
            volume.yOrigin = getYOrigin();

            if(zDomain != objectZDomain)
            {
                StsException.systemError("Error in StsProject.addToProject(StsSeismicVolume).\n" +
                        "Cannot add new volume to current boundingBoxes as volume domain is " +
                        StsParameters.TD_ALL_STRINGS[objectZDomain] +
                        " and project zDomain is " + StsParameters.TD_ALL_STRINGS[zDomain] + ".");
                return false;
            }
            rotatedBoundingBox.addVolumeBoundingBox(volume);
            unrotatedBoundingBox.addRotatedBoundingBox(rotatedBoundingBox);
            checkChangeZRange(volume, objectZDomain);
            adjustBoundingBoxes(true, false);
            //            runCompleteLoading();
            //setCursorDisplayXYAndGridCheckbox(false);  //how do you know this project is 3d? SWC 8/31/09
            setCursorDisplayXYAndGridCheckbox(false);  //reenable wrw
            resetCropVolume();
            //           currentModel.win3d.cursor3dPanel.setSliderValues();
            //            rangeChanged();
            objectPanelChanged();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.addToProject() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public void setCursorDisplayXYAndGridCheckbox(final boolean isXY)
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                StsWin3dFull[] windows = currentModel.getParentWindows();
                for(int n = 0; n < windows.length; n++)
                    windows[n].setCursorDisplayXYAndGridCheckbox(isXY);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    /**
     * This seismicBoundingBox bounds an unrotated seismic object such as a 2D line which has no inherent rotation angle
     * and need not be congruent with existing rotated objects (seismicVolumes and surfaces).
     * If an origin has already been set by some object, then compute x and y min and max relative to this origin to represent
     * an unrotated bounding box around this object.
     * If an origin has not been set, then set the origin from this unrotated seismic object.  Once fixed, the origin can't
     * be moved.
     * Subsequent to this operation if a rotatable seismic object is loaded it will use the same origin, but will set the
     * angle to something other than zero which redefines the local coordinate system.  Any unrotatedClasses will have
     * their appropriate coordinates recomputed in this new rotated coordinate system.  These unrotatedClasses will be
     * deleted from the current list indicating they have been dealt with.
     *
     * @param volume
     * @return
     */
    public boolean addUnrotatedBoxToProject(StsSeismic volume)
    {
        try
        {
            byte objectZDomain = StsParameters.getZDomainFromString(volume.zDomain);
            if(!checkSetZDomain(objectZDomain)) return false;

            unrotatedBoundingBox.addUnrotatedBoundingBox(volume);

            if(zDomain != objectZDomain)
            {
                StsException.systemError("Error in StsProject.addToProject(StsSeismicBoundingBox).\n" +
                        "Cannot add new volume to current boundingBoxes as volume domain is " +
                        StsParameters.TD_ALL_STRINGS[objectZDomain] +
                        " and project zDomain is " + StsParameters.TD_ALL_STRINGS[zDomain] + ".");
                return false;
            }

            checkChangeZRange(volume, objectZDomain);
            rangeChanged();
            objectPanelChanged();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.addToProject() failed.", e, StsException.WARNING);
            return false;
        }
    }

    /**
     * if we already have a volume loaded that we need to be congruent with in Z, adjust z range accordingly.
     * Otherwise, this volume will be considered the congruent basis for any subsequent volume load or velocity model construction.
     */
    private void checkChangeZRange(StsSeismicBoundingBox volume, byte zDomain)
    {
        checkChangeZRange(volume.zMin, volume.zMax, volume.zInc, zDomain);
        setIncLock(zDomain);
    }

    private void checkChangeZRange(float zMin, float zMax, float zInc, byte zDomain)
    {
        if(zInc == 0.0f)
        {
            StsException.systemError(this, "checkChangeZRange", "zInc cannot be zero.");
            return;
        }
        if(zDomain == TD_DEPTH)
        {
            if(depthIncLocked)
            {
                setDepthMin(StsMath.intervalRoundDown(zMin, depthMin, depthInc));
                setDepthMax(StsMath.intervalRoundUp(zMax, depthMin, depthInc));
            }
            else
            {
                float newDepthMin = Math.min(zMin, depthMin);
                setDepthInc(zInc);
                setDepthMin(StsMath.intervalRoundDown(newDepthMin, zMin, depthInc));
                float newDepthMax = Math.max(zMax, depthMax);
                setDepthMax(StsMath.intervalRoundUp(newDepthMax, zMin, depthInc));
                setDepthIncLock();
            }
        }
        else if(zDomain == TD_TIME)
        {
            if(timeIncLocked)
            {
                setTimeMin(StsMath.intervalRoundDown(zMin, timeMin, timeInc));
                setTimeMax(StsMath.intervalRoundUp(zMax, timeMin, timeInc));
            }
            else
            {
                float newTimeMin = Math.min(zMin, timeMin);
                setTimeInc(zInc);
                setTimeMin(StsMath.intervalRoundDown(newTimeMin, zMin, timeInc));
                float newTimeMax = Math.max(zMax, timeMax);
                setTimeMax(StsMath.intervalRoundUp(newTimeMax, zMin, timeInc));
                setTimeIncLock();
            }
        }
    }

    public void checkChangeZRange(float zMin, float zMax, byte zDomain)
    {
        if(zDomain == TD_DEPTH)
        {
            if(depthMin == StsParameters.largeFloat)
                setDepthMin(zMin);
            else
                setDepthMin(StsMath.intervalRoundDown(zMin, depthMin, depthInc));
            setDepthMax(StsMath.intervalRoundUp(zMax, depthMin, depthInc));

        }
        else if(zDomain == TD_TIME)
        {
            if(timeMin == StsParameters.largeFloat)
                setTimeMin(zMin);
            else
                setTimeMin(StsMath.intervalRoundDown(zMin, timeMin, timeInc));
            setTimeMax(StsMath.intervalRoundUp(zMax, timeMin, timeInc));
        }
    }

    public void checkChangeCongruentZMin(float zMin, byte zDomain)
    {
        if(zDomain == TD_DEPTH)
        {
            if(depthIncLocked)
            {
                setDepthMin(StsMath.intervalRoundDown(zMin, this.depthMin, this.depthInc));
            }
            else
            {
                float newDepthMin = Math.min(zMin, this.depthMin);
                setDepthMin(newDepthMin);
            }
        }
        else if(zDomain == TD_TIME)
        {
            if(timeIncLocked)
            {
                setTimeMin(StsMath.intervalRoundDown(zMin, this.timeMin, this.timeInc));
            }
            else
            {
                float newTimeMin = Math.min(zMin, this.timeMin);
                setTimeMin(newTimeMin);
            }
        }
    }

    public void checkChangeCongruentZMax(float zMax, byte zDomain)
    {
        checkChangeCongruentZMax(zMax, zDomain, false);
    }
    public void checkChangeCongruentZMax(float zMax, byte zDomain, boolean override)
    {
        if(zDomain == TD_DEPTH)
        {
            if(depthIncLocked)
            {
                setDepthMax(StsMath.intervalRoundUp(zMax, this.depthMin, this.depthInc), override);
            }
            else
            {
                float newDepthMax = Math.max(zMax, this.depthMax);
                setDepthMax(newDepthMax, override);
            }
        }
        else if(zDomain == TD_TIME)
        {
            if(timeIncLocked)
            {
                setTimeMin(StsMath.intervalRoundDown(zMax, this.timeMin, this.timeInc), override);
            }
            else
            {
                float newTimeMin = Math.min(zMax, this.timeMin);
                setTimeMin(newTimeMin, override);
                setTimeIncLock();
            }
        }
    }

    public void checkSetDepthRangeForVelocityModel()
    {
        if(velocityModel == null) return;
        velocityModel.checkSetDepthRange();
    }

    public void checkSetDepthRange(float maxAvgVelocity, float depthDatum, float timeDatum)
    {
        float newDepthMin = depthDatum + maxAvgVelocity * (timeMin - timeDatum);
        float newDepthMax = depthDatum + maxAvgVelocity * (timeMax - timeDatum);
        checkSetDepthRange(newDepthMin, newDepthMax);
    }

    public void checkSetDepthRange(float depthMin, float depthMax)
    {
        if(depthMax < this.depthMax) return;

        int approxNInc;
        if(timeInc != 0.0f)
            approxNInc = StsMath.ceiling((timeMax - timeMin) / timeInc);
        else
            approxNInc = 100;

        double[] scale = StsMath.niceScale(depthMin, depthMax, approxNInc, true);
        setDepthMin((float) scale[0]);
        setDepthMax((float) scale[1]);
        setDepthInc((float) scale[2]);
        objectPanelChanged();
    }

    public void checkSetDepthMax(float depthMax)
    {
        if(depthMax < this.depthMax) return;

        int approxNInc;
        if(timeInc != 0.0f)
            approxNInc = StsMath.ceiling((timeMax - timeMin) / timeInc);
        else
            approxNInc = 100;

        double[] scale = StsMath.niceScale(this.depthMin, depthMax, approxNInc, true);
        setDepthMin((float) scale[0]);
        setDepthMax((float) scale[1]);
        setDepthInc((float) scale[2]);
        objectPanelChanged();
    }

      private void setUserDepthMin(float depthMin)
      {
          if(this.depthMin == depthMin) return;
          if(depthMin >= depthMax)
          {
              StsException.systemDebug(this, "setDepthMin", "Problem:  depthMin( " + depthMin + " >= depthMax" + depthMax);
              return;
          }
          checkChangeCongruentZMin(depthMin, TD_DEPTH);
      }

      private void setDepthMin(float depthMin)
      {
          setDepthMin(depthMin, false);
      }

      private void setDepthMin(float depthMin, boolean override)
      {
          if(!override && this.depthMin <= depthMin) return;
          if(debug) StsException.systemDebug(this, "setDepthMin", "depthMin changed from " + this.depthMin + " to " + depthMin);
          this.depthMin = depthMin;
          if(zDomain == TD_DEPTH)
          {
              displayBoundingBox.setZMin(depthMin);
              rotatedBoundingBox.setZMin(depthMin);
              unrotatedBoundingBox.setZMin(depthMin);
              currentModel.cursorZRangeChanged();
          }
          this.dbFieldChanged("depthMin", depthMin);
          currentModel.win3dDisplay();
      }

      public float getDepthMin() { return depthMin; };

      public float getUserDepthMin() { return getDepthMin(); }

      private void setUserDepthMax(float depthMax)
      {
          if(this.depthMin == depthMin) return;
          if(depthMax <= depthMin)
          {
              StsException.systemDebug(this, "setDepthMin", "Problem:  depthMax " + depthMax + " <= depthMin" + depthMin);
              return;
          }
          checkChangeCongruentZMax(depthMax, TD_DEPTH);
      }

      private void setDepthMax(float depthMax)
      {
            setDepthMax(depthMax, false);
      }

      private void setDepthMax(float depthMax, boolean override)
      {
          if(!override && this.depthMax >= depthMax) return;
          if(debug) StsException.systemDebug(this, "setDepthMax", "depthMax changed from " + this.depthMax + " to " + depthMax);
          this.depthMax = depthMax;
          if(zDomain == TD_DEPTH)
          {
              displayBoundingBox.setZMax(depthMax);
              rotatedBoundingBox.setZMax(depthMax);
              unrotatedBoundingBox.setZMax(depthMax);
              currentModel.cursorZRangeChanged();
              setGridZ(depthMax);
          }

          this.dbFieldChanged("depthMax", depthMax);
          currentModel.win3dDisplay();
      }

      public float getDepthMax() { return depthMax; };

    public float getUserDepthMax() { return getDepthMax(); }

    /** this is an approximation used only for debug display purposes. */
    public double getRatioRangeZorT(double z, byte zDomainData, byte zDomainProject)
    {
        if(zDomainData == zDomainProject) return z;
        if(zDomainData == TD_DEPTH) // zDomainProject == TD_TIME
            return timeMin + (timeMax - timeMin)*(z - depthMin)/(depthMax - depthMin);
        else // zDomainData == TD_TIME && zDomainProject ==  TD_DEPTH
            return depthMin + (depthMax - depthMin)*(z - timeMin)/(timeMax - timeMin);
    }
    /*
      public void setDepthDatum(float datum)
      {
          if(datum == this.depthDatum) return;
          this.depthDatum = datum;
          this.dbFieldChanged("depthDatum", depthDatum);
      }
      public float getDepthDatum() { return depthDatum; }

      public void setTimeDatum(float datum)
      {
          if(datum == this.timeDatum) return;
          this.timeDatum = datum;
          this.dbFieldChanged("timeDatum", timeDatum);
      }
      public float getTimeDatum() { return timeDatum; }
       */
    public void resetProjectBounds()
    {
        adjustBoundingBoxes(true, true);
        rangeChanged();
        currentModel.win3d.cursor3d.initialize();
    }

    private void setTimeMin(float timeMin)
    {
        setTimeMin(timeMin, false);
    }

    private void setTimeMin(float timeMin, boolean override)
    {
        if(!override && timeMin >= this.timeMin) return;
        if(debug) StsException.systemDebug(this, "setTimeMin", "timeMin changed from " + this.timeMin + " to " + timeMin);
        this.timeMin = timeMin;
        if(zDomain == TD_DEPTH)
        {
              displayBoundingBox.setZMin(timeMin);
              rotatedBoundingBox.setZMin(timeMin);
              unrotatedBoundingBox.setZMin(timeMin);
              currentModel.cursorZRangeChanged();
        }
        this.dbFieldChanged("timeMin", timeMin);
    }

    public float getTimeMin() { return timeMin; }

    private void setTimeMax(float timeMax)
    {
        setTimeMax(timeMax, false);
    }

    private void setTimeMax(float timeMax, boolean override)
    {
        if(!override && timeMax <= this.timeMax) return;
        float oldTimeMax = this.timeMax;
        this.timeMax = timeMax;
        if(debug) StsException.systemDebug(this, "setTimeMax", "timeMax changed from " + oldTimeMax + " to " + this.timeMax);
        if(zDomain == TD_TIME)
        {
            displayBoundingBox.setZMax(timeMax);
            rotatedBoundingBox.setZMax(timeMax);
            unrotatedBoundingBox.setZMax(timeMax);
            currentModel.cursorZRangeChanged();
            setGridZ(timeMax);
        }
        this.dbFieldChanged("timeMax", timeMax);
    }

    public float getTimeMax() { return timeMax; }

    private void checkSetTimeInc(float newInc)
    {
        if(timeIncLocked || timeInc != 0.0f && newInc >= timeInc) return;

        if(zDomain == TD_TIME)
            rotatedBoundingBox.setZInc(newInc);
        if(timeInc == newInc) return;
        timeInc = newInc;

        dbFieldChanged("timeInc", timeInc);
    }

    private void checkSetDepthInc(float newInc)
    {
        if(depthIncLocked || depthInc != 1.0f && newInc >= depthInc)
            return; // hacque: default project zInc is 1.0; override.  Should have flag instead.

        if(timeInc != 0.0f)
        {
            int approxNInc = StsMath.ceiling((timeMax - timeMin) / timeInc);
            double[] scale = StsMath.niceScale(depthMin, depthMax, approxNInc, true);
            newInc = (float) scale[2];
        }
        setDepthInc(newInc);
    }

    private void setDepthInc(float depthInc)
    {
        int nDepthIncs = StsMath.ceiling((depthMax - depthMin) / depthInc);
        if(nDepthIncs > nZIncsMax)
        {
            double[] scale = StsMath.niceScale(depthMin, depthMax, nZIncsMax, true);
            depthInc = (float) scale[2];
        }
        depthInc = Math.max(depthInc, minZInc);
        if(zDomain == TD_DEPTH)
        {
            rotatedBoundingBox.zInc = depthInc;
        }
        if(this.depthInc == depthInc) return;
        this.depthInc = depthInc;
        dbFieldChanged("depthInc", depthInc);
    }

    private void setTimeInc(float timeInc)
    {
        timeInc = Math.max(timeInc, minZInc);
        if(zDomain == TD_TIME)
        {
            rotatedBoundingBox.zInc = timeInc;
        }
        if(this.timeInc == timeInc) return;
        this.timeInc = timeInc;
        dbFieldChanged("timeInc", timeInc);
    }

    public void setIncLock(byte zDomain)
    {
        if(zDomain == TD_DEPTH)
            setDepthIncLock();
        else if(zDomain == TD_TIME)
            setTimeIncLock();

    }
    public void setDepthIncLock()
    {
        if(depthIncLocked) return;
        fieldChanged("depthIncLocked", true);
    }

    public void setTimeIncLock()
    {
        if(timeIncLocked) return;
        fieldChanged("timeIncLocked", true);
    }


    public void rangeChanged()
    {
        resetCropVolume();
        if(currentModel.win3d.cursor3dPanel == null) return;
        currentModel.win3d.cursor3dPanel.setSliderValues();
    }

    public boolean addToProject(StsSurface surface)
    {
        byte objectZDomain = surface.getZDomainOriginal();
        if(!checkSetZDomain(objectZDomain))
        {
            return false;
        }

        boolean originChanged = checkSetOrigin(surface.xOrigin, surface.yOrigin);

        boolean angleAlreadySet = rotatedBoundingBox.getAngleSet();
        if(!checkSetAngle(surface.angle)) return false;

        // If angle not already set, then the current rotatedBoundingBox is same as the
        // project.unrotatedBoundingBox. So if we have set the angle now, we start with
        // a new unrotatedBoundingbox and add this rotated object to it.
        if(!angleAlreadySet)
        {
            rotatedBoundingBox.resetXYRange();
        }
        //        StsRotatedGridBoundingBox boundingBox = getRotatedBoundingBox();

        float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY(surface.xOrigin, surface.yOrigin);
        surface.xMin = xy[0];
        surface.yMin = xy[1];
        surface.xMax = surface.xMin + (surface.nCols - 1) * surface.xInc;
        surface.yMax = surface.yMin + (surface.nRows - 1) * surface.yInc;
        surface.xOrigin = getXOrigin();
        surface.yOrigin = getYOrigin();

        if(zDomain != objectZDomain)
        {
            StsException.systemError("Error in StsProject.addToProject(StsSurface).\n" +
                    "Cannot add new surface to current boundingBoxes as surface domain is " +
                    StsParameters.TD_ALL_STRINGS[objectZDomain] +
                    " and project zDomain is " + StsParameters.TD_ALL_STRINGS[zDomain] + ".");
            return false;
        }
        rotatedBoundingBox.addVolumeBoundingBox(surface);
        rotatedBoundingBox.checkMakeCongruent(surface);
        //        rotatedBoundingBox.adjustRowColNumbering(surface);
        float zDif = surface.zMax - surface.zMin;
        double[] scale = niceScale(surface.zMin - zDif / 2, surface.zMax + zDif / 2);
        float newMin = (float) scale[0];
        float newMax = (float) scale[1];
        float newInc = (float) scale[2];
        boolean changed = false;
        if(zDomain == TD_DEPTH)
        {
            if(newMin < depthMin)
            {
                changed = true;
                setDepthMin(newMin);
            }
            if(newMax > depthMax)
            {
                changed = true;
                setDepthMax(newMax);
            }
            if(changed) checkSetDepthInc(newInc);
        }
        else if(zDomain == TD_TIME)
        {
            if(newMin < timeMin)
            {
                changed = true;
                setTimeMin(newMin);
            }
            if(newMax > timeMax)
            {
                changed = true;
                setTimeMax(newMax);
            }
            if(changed) checkSetTimeInc(newInc);
        }
        //        if(changed)
        {
            //        rotatedBoundingBox.setZRange((float)scale[0], (float)scale[1], (float)scale[2]);
            //        rotatedBoundingBox.adjustZRange(surface.zMin - zDif / 2, surface.zMax + zDif / 2);
            //        boundingBox.addVolumeBoundingBox(this);
            StsBoundingBox unrotatedBoundingBox = getUnrotatedBoundingBox();
            unrotatedBoundingBox.addRotatedBoundingBox(rotatedBoundingBox);
            rangeChanged();
        }
        surface.setRelativeRotationAngle();
        objectPanelChanged();
        return true;
    }

    /**
     * Get the Project minimum in X
     *
     * @return minimum X
     */
    public float getXMin()
    {
        return rotatedBoundingBox.xMin;
    }

    /**
     * Get Project maximum in X
     *
     * @return maximum X
     */
    public float getXMax()
    {
        return rotatedBoundingBox.xMax;
    }

    /**
     * Get the Project size in X dimension
     *
     * @return size in X
     */
    public float getXSize()
    {
        return rotatedBoundingBox.xMax - rotatedBoundingBox.xMin;
    }

    /**
     * Get the X increment - actual
     *
     * @return X increment
     */
    public float getXInc()
    {
        return rotatedBoundingBox.xInc;
    }

    /**
     * Get the project grid displayed x increment (dotted line spacing in project grid)
     *
     * @return project grid displayed x increment
     */
    public float getGridDX()
    {
        return displayBoundingBox.gridDX;
    }

    /**
     * Get absolute X origin
     *
     * @return absolute X origin
     */
    public double getXOrigin()
    {
        return rotatedBoundingBox.xOrigin;
    }

    /**
     * Get X at center of Project
     *
     * @return center of gravity X
     */
    public float getXCenter()
    {
        return rotatedBoundingBox.getXCenter();
    }

    /**
     * Get the Project minimum in Y
     *
     * @return minimum y
     */
    public float getYMin()
    {
        return rotatedBoundingBox.yMin;
    }

    /**
     * Get Project maximum in Y
     *
     * @return maximum Y
     */
    public float getYMax()
    {
        return rotatedBoundingBox.yMax;
    }

    /** get slice coordinate from local z coordinate */
    public int getNearestSliceCoor(boolean isDepth, float z)
    {
        int slice;
        if(isDepth)
            slice = Math.round((z - depthMin) / depthInc);
        else
            slice = Math.round((z - timeMin) / timeInc);
        if(slice < 0 || slice >= rotatedBoundingBox.nSlices) return -1;
        return slice;
    }

    public String toString()
    {
        return "Project";
    }

    public boolean isPlanar()
    {
        return true;
    }

	public void adjustRotatedBoundingBoxGrid(StsRotatedGridBoundingBox inputBoundingBox, byte zDomain)
	{
		// hack to fix rotatedBoundingBox
		rotatedBoundingBox.setRowColIndexRanges();
		// define boundingBox xInc and yInc so that they are integral values of current xInc or yInc or vice versa
		int rowMin = rotatedBoundingBox.getFloorRowCoor(inputBoundingBox.yMin);
		inputBoundingBox.yMin = rotatedBoundingBox.getYCoor(rowMin);
		int rowMax = rotatedBoundingBox.getCeilingRowCoor(inputBoundingBox.yMax);
		inputBoundingBox.yMax = rotatedBoundingBox.getYCoor(rowMax);
		inputBoundingBox.yInc = rotatedBoundingBox.yInc;
		inputBoundingBox.nRows = rowMax - rowMin + 1;

		int colMin = rotatedBoundingBox.getFloorColCoor(inputBoundingBox.xMin);
		inputBoundingBox.xMin = rotatedBoundingBox.getXCoor(colMin);
		int colMax = rotatedBoundingBox.getCeilingColCoor(inputBoundingBox.xMax);
		inputBoundingBox.xMax = rotatedBoundingBox.getXCoor(colMax);
		inputBoundingBox.xInc = rotatedBoundingBox.xInc;
		inputBoundingBox.nCols = colMax - colMin + 1;
		
		if(zDomain == TD_TIME && supportsTime())
            inputBoundingBox.resetZRange(timeMin, timeMax, timeInc);
        else if(zDomain == TD_DEPTH && supportsDepth())
            inputBoundingBox.resetZRange(depthMin, depthMax, depthInc);
	}

    public void setCellSize(float size)
    {
    	this.cellSize = size;
    }

    public float getCellSize()
    {
    	return cellSize;
    }

    /** compute a gridSize based on the gridSize provided such that it is an integral fraction of the current gridIncrement
     *  or an integral multiplier of the current gridIncrement subject to the limit of a maxGridSize
     * @param inc the current grid increment
     * @param gridSize desired grid cell size
     * @param gridLength length of the current grid in this direction
     * @return resulting grid cell size
     */
  	static public float computeIntegralGridSize(float inc, float gridSize, float gridLength)
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
    /**
     * Get the Project Grid Frequency as String
     *
     * @return grid frequency in model units
     */
    public String getGridFrequencyString()
    {
        if(gridFrequency == AUTO)
        {
            return new String("Automatic");
        }
        else
        {
            return new Float(gridFrequency).toString();
        }
    }

    /**
     * Get the Project Grid Frequency as String
     *
     * @return grid frequency in model units
     */
    public String getZGridFrequencyString()
    {
        if(zGridFrequency == AUTO)
        {
            return new String("Automatic");
        }
        else
        {
            return new Float(zGridFrequency).toString();
        }
    }

    /**
     * Get Project Label Frequency String
     *
     * @return label frequency in model units
     */
    public String getLabelFrequencyString()
    {
        if(labelFrequency == AUTO)
        {
            return new String("Automatic");
        }
        else
        {
            return new Float(labelFrequency).toString();
        }
    }

    /**
     * Get the Project Grid Frequency
     *
     * @return grid frequency in model units
     */
    public float getGridFrequency()
    {
        return gridFrequency;
    }

    /**
     * Get Project Label Frequency
     *
     * @return label frequency in model units
     */
    public float getLabelFrequency()
    {
        return labelFrequency;
    }

    /** Set the Project Z Grid Frequency */
    public void setZGridFrequency(float freq)
    {
        zGridFrequency = freq;
        // dbFieldChanged("zGridFrequency", freq);
        return;
    }

    /** Set the Project Grid Frequency */
    public void setGridFrequency(float freq)
    {
        gridFrequency = freq;
        // dbFieldChanged("gridFrequency", freq);
        return;
    }

    /** Set Project Label Frequency */
    public void setLabelFrequency(float freq)
    {
        labelFrequency = freq;
        //		setDisplayField("labelFrequency", freq);
        //		dbFieldChanged("labelFrequency", freq);
        return;
    }

    /** Set the Project Grid Frequency String */
    public void setGridFrequencyString(String freq)
    {
        if(freq.equals("Auto") || freq.equals("Automatic"))
        {
            gridFrequency = AUTO;
        }
        else
        {
            gridFrequency = (new Float(freq)).floatValue();
            if((labelFrequency % gridFrequency != 0) && (labelFrequency != AUTO))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "Label frequency must be a multiple of grid frequency.");
                setGridFrequency(AUTO);
            }
            currentModel.win3dDisplayAll();
        }
        // dbFieldChanged("gridFrequency", gridFrequency);
        return;
    }

    /** Set the Project Grid Frequency String */
    public void setZGridFrequencyString(String freq)
    {
        if(freq.equals("Auto") || freq.equals("Automatic"))
        {
            zGridFrequency = AUTO;
        }
        else
        {
            zGridFrequency = (new Float(freq)).floatValue();
            if((labelFrequency % zGridFrequency != 0) && (labelFrequency != AUTO))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "Label frequency must be a multiple of grid frequency.");
                setGridFrequency(AUTO);
            }
            currentModel.win3dDisplayAll();
        }
        // dbFieldChanged("zGridFrequency", zGridFrequency);
        return;
    }

    public float getZGridFrequency()
    {
        return zGridFrequency;
    }

    /** Set Project Label Frequency String */
    public void setLabelFrequencyString(String freq)
    {
        if(freq.equals("Auto") || freq.equals("Automatic"))
        {
            labelFrequency = AUTO;
        }
        else
        {
            labelFrequency = (new Float(freq)).floatValue();
            if((labelFrequency % gridFrequency != 0) && (gridFrequency != AUTO))
            {
                new StsMessage(currentModel.win3d, StsMessage.WARNING, "Label frequency must be a multiple of grid frequency.");
                setLabelFrequency(AUTO);
            }
            currentModel.win3dDisplayAll();
        }
        // dbFieldChanged("labelFrequency", labelFrequency);
        return;
    }

    /**
     * Get the Project size in Y dimension
     *
     * @return size in Y
     */
    public float getYSize()
    {
        return rotatedBoundingBox.yMax - rotatedBoundingBox.yMin;
    }

    public float getSize()
    {
        return Math.max(getXSize(), getYSize());
    }

    /**
     * Get the Y increment - actual
     *
     * @return Y increment
     */
    public float getYInc()
    {
        return rotatedBoundingBox.yInc;
    }

    /**
     * Get the project grid displayed y increment (dotted line spacing in project grid)
     *
     * @return project grid displayed y increment
     */
    public float getGridDY()
    {
        return displayBoundingBox.gridDY;
    }

    /**
     * Get absolute Y origin
     *
     * @return absolute Y origin
     */
    public double getYOrigin()
    {
        return rotatedBoundingBox.yOrigin;
    }

    /**
     * Get Y at center of Project
     *
     * @return center of gravity Y
     */
    public float getYCenter()
    {
        return rotatedBoundingBox.getYCenter();
    }

    //    public float getAngle() { return angle; }

    /**
     * Get the Project minimum in Z
     *
     * @return minimum Z
     */
    public float getZorTMin()
    {
        if(isDepth)
            return depthMin;
        else
            return timeMin;
    }

    public float getZorTMin(boolean isDepth)
    {
        if(isDepth)
            return depthMin;
        else
            return timeMin;
    }

    /**
     * Get Project maximum in Z
     *
     * @return maximum Z
     */
    public float getZorTMax()
    {
        if(isDepth)
        {
            return depthMax;
        }
        else
        {
            return timeMax;
        }
    }

    public float getZorTMax(boolean isDepth)
    {
        if(isDepth)
            return depthMax;
        else
            return timeMax;
    }

    /**
     * Get the Z increment - actual
     *
     * @return Z increment
     */
    public float getZorTInc()
    {
        if(isDepth)
            return depthInc;
        else
            return timeInc;
    }

    public float getZorTInc(boolean isDepth)
    {
        if(isDepth)
            return depthInc;
        else
            return timeInc;
    }

    public float getGridZ()
    {
        if(displayBoundingBox.gridZ != StsParameters.nullValue) return displayBoundingBox.gridZ;
        if(isDepth)
            return depthMax;
        else
            return timeMax;
    }

    public float getCultureZ(float x, float y)
    {
        return getGridZ();
    }

    /**
     * Get angle from global +X to project +X axes
     *
     * @return angle
     */
    public float getAngle()
    {
        return rotatedBoundingBox.angle;
    }

    public float getDisplayXMin()
    {
        return displayBoundingBox.xMin;
    }

    public float getDisplayXMax()
    {
        return displayBoundingBox.xMax;
    }

    public float getDisplayYMin()
    {
        return displayBoundingBox.yMin;
    }

    public float getDisplayYMax()
    {
        return displayBoundingBox.yMax;
    }
/*
    public float getDepthMin()
    {
        if(supportsDepth())
            return depthMin;
        else
            return 0.0f;
    }

    public float getDepthMax()
    {
        if(supportsDepth())
            return depthMax;
        else
            return 0.0f;
    }
*/
    public float getDepthInc()
    {
        return depthInc;
    }

    public long getProjectTime()
    {
        return projectTime;
    }

    public String getProjectTimeAsString()
    {
        Date date = new Date(projectTime);
        SimpleDateFormat format = new SimpleDateFormat(getTimeDateFormatString());
        String time = format.format(date);
        return time;
    }

    public boolean setProjectTime(long time, boolean updateTb)
    {
        setProjectTime(time);
        if(updateTb)
        {
            StsWin3d window = currentModel.win3d;
            StsTimeActionToolbar tb = (StsTimeActionToolbar) window.getTimeActionToolbar();
            tb.stopTime();
            if(tb != null)
                tb.updateTime();
        }
        return true;
    }

    public void stopProjectTime()
    {
        StsWin3d window = currentModel.win3d;
        StsTimeActionToolbar tb = (StsTimeActionToolbar) window.getTimeActionToolbar();
        if(tb == null) return;
        if(tb.isRunning())
            tb.stopTime();
    }

    public boolean setProjectTime(long time)
    {
        boolean success = true;
        projectTime = time;
        //if(projectTime >= System.currentTimeMillis())
        //    success = false;
        //else
        //    success = true;

        // Poll all monitor objects. The polling frequency is based on the tb clock update rate.
        // Even if not in real-time mode, if the clock is running monitor objects are being polled.
        // This way the data is added to the monitored objects but the time is still applied when displayed
        // so the new data may not be visible yet.
        StsObject[] monitorObjects = currentModel.getObjectList(StsMonitor.class);
        for(int i = 0; i < monitorObjects.length; i++)
            ((StsMonitor) monitorObjects[i]).poll();

        currentModel.viewObjectRepaint(this, this);
        return success;
    }

    public boolean incrementProjectTime(long mseconds)
    {
        return setProjectTime(projectTime + mseconds);
    }

    public void setProjectTimeToCurrentTime(boolean updateToolbar)
    {
        setProjectTime(System.currentTimeMillis(), updateToolbar);
    }

    public void setProjectTimeToCurrentTime()
    {
        projectTime = System.currentTimeMillis();
        currentModel.viewObjectRepaint(this, this);
    }

    public void clearProjectTimeDuration()
    {
        projectTimeDuration = 0;
        currentModel.viewObjectRepaint(this, this);
        return;
    }

    public long getProjectTimeDuration()
    {
        return projectTimeDuration;
    }

    public String getTimeString(long time)
    {
        return getTimeDateFormat().format(new Date(time));
    }

    public long getTime(String timeString)
    {
        long nTime = -1l;
        try
        {
            nTime = getTimeDateFormat().parse(timeString).getTime();
        }
        catch(Exception ex)
        {
            System.out.println("StsProject:Failed to parse supplied time string.");
        }
        return nTime;
    }

    public boolean setProjectTimeDuration(long startTime, long endTime)
    {
        projectTimeDuration = endTime - startTime;
        StsWin3d window = currentModel.win3d;
        StsTimeActionToolbar tb = (StsTimeActionToolbar) window.getTimeActionToolbar();
        if(tb == null) return true;
        if(tb.isRunning())
            tb.stopTime();
        return true;
    }

    public boolean isTimeRunning()
    {
        StsWin3d window = currentModel.win3d;
        StsTimeActionToolbar tb = (StsTimeActionToolbar) window.getTimeActionToolbar();
        return tb.isRunning();
    }

    public boolean isRealtime()
    {
        StsWin3d window = currentModel.win3d;
        if(window == null)
            return false;
        StsTimeActionToolbar tb = (StsTimeActionToolbar) window.getTimeActionToolbar();
        if(tb != null)
            return tb.isRealtime();
        else
            return false;
    }

 /*
    public float getTimeMin()
    {
        if(supportsTime())
            return timeMin;
        else
            return 0.0f;
    }

    public float getTimeMax()
    {
        if(supportsTime())
            return timeMax;
        else
            return 0.0f;
    }
*/
    public float getTimeInc()
    {
        return timeInc;
    }

    /**
     * Set the background color of the graphics environment
     *
     * @argument Java Color
     */
    public void setBackgroundColor(StsColor color)
    {
        if(backgroundColor == color) return;
        backgroundColor = new StsColor(color);
        //        currentModel.properties.set(backgroundColorBean.getBeanKey(), backgroundColor.toString());
        currentModel.win3dDisplayAll();
    }

    /**
     * Get the current background color of the graphics environment
     *
     * @return Java Color
     */
    public StsColor getBackgroundColor()
    {
        return backgroundColor;
    }

    /**
     * get current foreground color, which is inverse of background color
     *
     * @return StsColor
     */
    public StsColor getForegroundColor()
    {
        return StsColor.getInverseStsColor(backgroundColor);
    }

    /**
     * Get the current background color of the graphics environment as StsColor
     *
     * @return StsColor
     */
    public StsColor getBackgroundStsColor()
    {
        return backgroundColor;
    }

    /**
     * Set the grid color of the graphics environment
     *
     * @argument Java Color
     */
    public void setGridColor(StsColor color)
    {
        if(gridColor == color) return;
        gridColor = new StsColor(color);
        //		currentModel.properties.set(gridColorBean.getBeanKey(), gridColor.toString());
        currentModel.win3dDisplayAll();
    }

    /**
     * Get the current grid color of the graphics environment
     *
     * @return Java Color
     */
    public StsColor getGridColor()
    {
        return gridColor;
    }

    /**
     * Get the current grid color of the graphics environment as StsColor
     *
     * @return StsColor
     */
    public StsColor getStsGridColor()
    {
        return gridColor;
    }

    /**
     * Set the color of the COG in graphics environment
     *
     * @argument Java Color
     */
    public void setCogColor(StsColor color)
    {
        if(cogColor == color) return;
        cogColor = new StsColor(color);
        //		currentModel.properties.set(cogColorBean.getBeanKey(), cogColor.toString());
        currentModel.win3dDisplayAll();
    }

    /**
     * Set the color of the timing lines in graphics environment
     *
     * @argument Java Color
     */
    public void setTimingColor(StsColor color)
    {
        if(timingColor == color) return;
        timingColor = new StsColor(color);
        //		currentModel.properties.set(timingColorBean.getBeanKey(), timingColor.toString());
        currentModel.win3dDisplayAll();
    }
    /**
     * Set the color of the timing lines in graphics environment
     *
     * @argument Java Color
     */
    /*
        public void setWellPanelColor(StsColor color)
		{
			if(wellPanelColor == color) return;
			wellPanelColor = color;
	//		currentModel.properties.set(timingColorBean.getBeanKey(), timingColor.toString());
			currentModel.win3dDisplayAll();
		}
      */

    /**
     * Get the current COG color of the graphics environment
     *
     * @return Java Color
     */
    public StsColor getCogColor()
    {
        return cogColor;
    }

    /**
     * Get the current COG color of the graphics environment as StsColor
     *
     * @return StsColor
     */
    public StsColor getStsCogColor()
    {
        return cogColor;
    }

    /**
     * Get the current timing line color of the graphics environment
     *
     * @return Java Color
     */
    public StsColor getTimingColor()
    {
        return timingColor;
    }
    /**
     * Get the current timing line color of the graphics environment
     *
     * @return Java Color
     */
    /*
        public StsColor getWellPanelColor()
		{
			return wellPanelColor;
		}

		public StsColor getStsWellPanelColor()
		{
			return wellPanelColor;
		}
       */

    /**
     * Get the current grid color of the graphics environment as StsColor
     *
     * @return StsColor
     */
    public StsColor getStsTimingColor()
    {
        return timingColor;
    }

    /** Set the minimum Z
     * @params minimum z value
     */
    /*
      public void setZMin(float z)
      {
          rotatedBoundingBox.zMin = z;
      }
   */
    /** Set the maximum Z
     * @params maximum z value
     */
    /*
      public void setZMax(float z)
      {
          rotatedBoundingBox.zMax = z;
      }
  */
    /** Set the Z increment
     * @params z increment
     * If this is actually called, there may be confusion if both time and
     * depth are supported..check.
     */
    /*
          public void setZInc(float dz)
          {
              rotatedBoundingBox.zInc = dz;

              if (zDomain == TD_DEPTH)
               {
                   depthInc = dz;
               }
               else
               {
                   timeInc = dz;
              }
          }
      */

    /**
     * Set the Z grid increment
     *
     * @params z grid increment
     */

    public void setGridZ(float gridZ)
    {
        displayBoundingBox.gridZ = gridZ;
    }

    /**
     * Show the intersections between the cursor slices
     *
     * @params true to show
     */
    public void setShowIntersection(boolean showIntersection)
    {
        if(this.showIntersection == showIntersection) return;
        this.showIntersection = showIntersection;
        //		setDisplayField("showIntersection", showIntersection);
        //		dbFieldChanged("showIntersection", showIntersection);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine whether the intesections are to be shown
     *
     * @returns true for show
     */
    public boolean getShowIntersection()
    {
        return showIntersection;
    }

    /**
     * Set the Project to visible
     *
     * @params true to show
     */
    public void setIsVisible(boolean isVisible)
    {
        if(this.isVisible == isVisible) return;
        this.isVisible = isVisible;
        //		setDisplayField("isVisible", isVisible);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine if the Project is visible
     *
     * @returns true for visible
     */
    public boolean getIsVisible()
    {
        return isVisible;
    }

    /**
     * Make the grid visible
     *
     * @params true to show
     */
    public void setShowGrid(boolean showGrid)
    {
        if(this.showGrid == showGrid) return;
        this.showGrid = showGrid;
        //		setDisplayField("showGrid", showGrid);
        //		dbFieldChanged("showGrid", showGrid);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine if the Grid is visible
     *
     * @returns true for visible
     */
    public boolean getShowGrid()
    {
        return showGrid;
    }

    /**
     * Make the grid on 3D views visible
     *
     * @params true to show
     */
    public void setShow3dGrid(boolean show3dGrid)
    {
        if(this.show3dGrid == show3dGrid) return;
        this.show3dGrid = show3dGrid;
        //		setDisplayField("show3dGrid", show3dGrid);
        dbFieldChanged("show3dGrid", show3dGrid);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine if the grid on 3D views are visible
     *
     * @returns true for visible
     */
    public boolean getShow3dGrid()
    {
        return show3dGrid;
    }

    /**
     * Make the grid on 2D views visible
     *
     * @params true to show
     */
    public void setShow2dGrid(boolean show2dGrid)
    {
        if(this.show2dGrid == show2dGrid) return;
        this.show2dGrid = show2dGrid;
        //		setDisplayField("show2dGrid", show3dGrid);
        dbFieldChanged("show2dGrid", show2dGrid);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine if the grid on 2D views are visible
     *
     * @returns true for visible
     */
    public boolean getShow2dGrid()
    {
        return show2dGrid;
    }

    /**
     * Make the labels visible
     *
     * @params true to show
     */
    public void setShowLabels(boolean showLabels)
    {
        if(this.showLabels == showLabels) return;
        this.showLabels = showLabels;
        //		setDisplayField("showLabels", showLabels);
        //		dbFieldChanged("showLabels", showLabels);
        currentModel.win3dDisplayAll();
    }

    /**
     * Determine if the labels are visible
     *
     * @returns true for visible
     */
    public boolean getShowLabels()
    {
        return showLabels;
    }

    /**
     * Set the value used to represent a null in map data
     *
     * @params null value
     */
    public void setMapGenericNull(float nullValue)
    {
        mapGenericNull = nullValue;
    }

    /**
     * Get the value used to represent a null in map data
     *
     * @returns map null value
     */
    public float getMapGenericNull()
    {
        return mapGenericNull;
    }

    /**
     * Set the value used to represent a null in log data
     *
     * @params null value
     */
    public void setLogNull(float logNull)
    {
        this.logNull = logNull;
    }

    /**
     * Get the value that represents a null in log data
     *
     * @returns log null value
     */
    public float getLogNull()
    {
        return logNull;
    }

    /** Set the min value used for time and depth increment for the project; i.e., they can't be smaller than this number */
    public void setMinZInc(float min)
    {
        if(minZInc == min) return;
        minZInc = min;
        dbFieldChanged("minZInc", minZInc);
    }

    /** Get the value that represents a null in log data */
    public float getMinZInc()
    {
        return minZInc;
    }

    /**
     * Get the Z domain
     *
     * @returns the Z domain item (defaults to time)
     */
    public byte getZDomain()
    {
        return zDomain;
    }

    public boolean isModelZDomain()
    {
        return modelZDomain == TD_NONE || modelZDomain == zDomain;
    }

    public boolean setToModelZDomain()
    {
        if(isModelZDomain()) return false;
        setZDomainAndRebuild(modelZDomain);
        return true;
    }

    public String getTimeUnitString()
    {
        return StsParameters.getTimeUnitString(timeUnits);
    }

    public void setTimeUnits(byte units)
    {
        timeUnits = units;
        dbFieldChanged("timeUnits", timeUnits);
    }

    public byte getTimeUnits()
    {
        return timeUnits;
    }

    public float getTimeScalar(byte from)
    {
        if(from == StsParameters.TIME_NONE)
            return 1.0f;
        else if(from == StsParameters.DIST_NONE)
            return 1.0f;
        else
            return StsParameters.TIME_SCALES[timeUnits] / StsParameters.TIME_SCALES[from];
    }

    public float getTimeScalar(String fromString)
    {
        byte from = StsParameters.getTimeUnitsFromString(fromString);
        return getTimeScalar(from);
    }

    public float calculateVelScaleMultiplier(String velocityUnits)
    {
        float scaleMultiplier = 1.0f;
        if(velocityUnits == StsParameters.VEL_UNITS_NONE)
        {
            return scaleMultiplier;
        }
        String projectUnits = this.getVelocityUnits();
        if(projectUnits == velocityUnits)
        {
            scaleMultiplier = 1.0f;
        }
        else if(projectUnits == StsParameters.VEL_M_PER_MSEC)
        {
            if(velocityUnits == StsParameters.VEL_M_PER_SEC)
                scaleMultiplier = 0.001f;
            else if(velocityUnits == StsParameters.VEL_FT_PER_MSEC)
                scaleMultiplier = 1.0f / StsParameters.DIST_FEET_SCALE;
            else if(velocityUnits == StsParameters.VEL_FT_PER_SEC)
                scaleMultiplier = 0.001f / StsParameters.DIST_FEET_SCALE;
        }
        else if(projectUnits == StsParameters.VEL_FT_PER_SEC)
        {
            if(velocityUnits == StsParameters.VEL_FT_PER_MSEC)
                scaleMultiplier = 1000;
            else if(velocityUnits == StsParameters.VEL_M_PER_MSEC)
                scaleMultiplier = 1000 * StsParameters.DIST_FEET_SCALE;
            else if(velocityUnits == StsParameters.VEL_M_PER_SEC)
                scaleMultiplier = StsParameters.DIST_FEET_SCALE;
        }
        else if(projectUnits == StsParameters.VEL_FT_PER_MSEC)
        {
            if(velocityUnits == StsParameters.VEL_FT_PER_SEC)
                scaleMultiplier = 0.001f;
            else if(velocityUnits == StsParameters.VEL_M_PER_MSEC)
                scaleMultiplier = StsParameters.DIST_FEET_SCALE;
            else if(velocityUnits == StsParameters.VEL_M_PER_SEC)
                scaleMultiplier = StsParameters.DIST_FEET_SCALE / 1000;
        }
        return scaleMultiplier;
    }

    public String getDepthUnitString()
    {
        return StsParameters.getDistanceUnitString(depthUnits);
    }

    public void setDepthUnits(byte units)
    {
        depthUnits = units;
        dbFieldChanged("depthUnits", depthUnits);
    }

    public byte getDepthUnits()
    {
        return depthUnits;
    }

    public float getDepthScalar(byte from)
    {
        if(from == StsParameters.DIST_NONE)
            return 1.0f;
        else
            return StsParameters.DIST_SCALES[depthUnits] / StsParameters.DIST_SCALES[from];
    }

    public float getDepthScalar(String fromString)
    {
        byte from = StsParameters.getDistanceUnitsFromString(fromString);
        return getDepthScalar(from);
    }

    public String getXyUnitString()
    {
        return StsParameters.getDistanceUnitString(xyUnits);
    }

    public void setXyUnits(byte units)
    {
        xyUnits = units;
        dbFieldChanged("xyUnits", xyUnits);
    }

    public byte getXyUnits()
    {
        return xyUnits;
    }

    public float getXyScalar(byte from)
    {
        if(from == StsParameters.DIST_NONE)
            return 1.0f;
        else
            return StsParameters.DIST_SCALES[xyUnits] / StsParameters.DIST_SCALES[from];
    }

    public float getXyScalar(String fromString)
    {
        byte from = StsParameters.getDistanceUnitsFromString(fromString);
        return getXyScalar(from);
    }

    public String getVerticalUnitsString()
    {
        if(zDomain == TD_DEPTH) return StsParameters.DIST_STRINGS[depthUnits];
        else                    return StsParameters.TIME_STRINGS[timeUnits];
    }

    public String getVerticalUnitsString(byte zDomain)
    {
        if(zDomain == TD_DEPTH) return StsParameters.DIST_STRINGS[depthUnits];
        else                    return StsParameters.TIME_STRINGS[timeUnits];
    }

    public void setTimeSeriesBkgdColor(StsColor color)
    {
        if(timeSeriesBkgdColor.equals(color)) return;
        timeSeriesBkgdColor = color;
        currentModel.viewObjectChanged(this, this);
    }

    public StsColor getTimeSeriesBkgdColor()
    {
        return timeSeriesBkgdColor;
    }

    public void setTimeSeriesFrgdColor(StsColor color)
    {
        if(timeSeriesFrgdColor.equals(color)) return;
        timeSeriesFrgdColor = color;
        currentModel.viewObjectChanged(this, this);
    }

    public StsColor getTimeSeriesFrgdColor()
    {
        return timeSeriesFrgdColor;
    }

    public void setTimeSeriesGridColor(StsColor color)
    {
        if(timeSeriesGridColor.equals(color)) return;
        timeSeriesGridColor = color;
        currentModel.viewObjectChanged(this, this);
    }

    public StsColor getTimeSeriesGridColor()
    {
        return timeSeriesGridColor;
    }

    public void setEnableTsGrid(boolean enable)
    {
        if(this.showTimeSeriesGrid == enable) return;
        this.showTimeSeriesGrid = enable;
        currentModel.viewObjectChanged(this, this);
    }

    public boolean getEnableTsGrid()
    {
        return showTimeSeriesGrid;
    }
    public void setCursorPositionByMouse(boolean enable)
    {
        if(this.cursorPositionByMouse == enable) return;
        this.cursorPositionByMouse = enable;
        currentModel.viewObjectChanged(this, this);
    }

    public boolean getCursorPositionByMouse()
    {
        return cursorPositionByMouse;
    }
    public String getVelocityUnits()
    {
        if(depthUnits == StsParameters.DIST_METER)
        {
            if(timeUnits == StsParameters.TIME_MSECOND)
                return StsParameters.VEL_M_PER_MSEC;
            else if(timeUnits == StsParameters.TIME_SECOND)
                return StsParameters.VEL_M_PER_SEC;
            else
                return StsParameters.VEL_UNITS_NONE;
        }
        else if(depthUnits == StsParameters.DIST_FEET)
        {
            if(timeUnits == StsParameters.TIME_MSECOND)
                return StsParameters.VEL_FT_PER_MSEC;
            else if(timeUnits == StsParameters.TIME_SECOND)
                return StsParameters.VEL_FT_PER_SEC;
            else
                return StsParameters.VEL_UNITS_NONE;
        }
        else
            return StsParameters.VEL_UNITS_NONE;
    }

    public void setDefaultSensorColor(Color color)
    {
        setDefaultSensorStsColor(new StsColor(color));
    }

    public Color getDefaultSensorColor()
    {
        return defaultSensorColor.getColor();
    }

    public void setDefaultSensorStsColor(StsColor color)
    {
        defaultSensorColor = new StsColor(color);
    }

    public StsColor getDefaultSensorStsColor()
    {
        return defaultSensorColor;
    }

    public void setDefaultWellColor(Color color)
    {
        setDefaultWellStsColor(new StsColor(color));
    }

    public StsColor getDefaultWellColor()
    {
        return defaultWellColor;
    }

    public void setDefaultWellStsColor(StsColor color)
    {
        defaultWellColor = new StsColor(color);
    }

    public StsColor getDefaultWellStsColor()
    {
        return defaultWellColor;
    }

    /** should be called only when we are toggling between time and depth states */
    public void setZDomainString(String zDomainString)
    {
        setZDomainAndRebuild(StsParameters.getZDomainFromString(zDomainString));
        //        resetProjectView();
    }

    public void setZDomain()
    {
        if(zDomainSupported != TD_TIME_DEPTH)
        {
            new StsMessage(currentModel.win3d, StsMessage.WARNING,
                    "Unable to switch domain without velocity information");
        }
        else
        {
            if(zDomain == TD_DEPTH)
            {
                setZDomainAndRebuild(TD_TIME);
            }
            else
            {
                setZDomainAndRebuild(TD_DEPTH);
            }
        }

    }

    protected void initializeZDomainAndRebuild()
    {
        setIsDepth(zDomain == TD_DEPTH);

        //		convertGridZ(oldZDomain);

        setBoundingBoxesZRanges();
        resetCropVolume();
        if(!windowInitialized()) return;
        //This looks redundant (z slider adjusted above and resetCropVolume() called).  Tom 12/12/07
        rangeChanged();
        currentModel.resetAllSliderValues();
        if(!isModelBuilt()) return;
        rebuildModel();
    }

    public boolean isModelBuilt()
    {
        StsBuiltModel builtModel = (StsBuiltModel) currentModel.getCurrentObject(StsBuiltModel.class);
        return builtModel != null;
    }

    private void rebuildModel()
    {
        final StsStatusPanel statusPanel = StsStatusPanel.constructStatusDialog(currentModel.win3d, "Rebuilding model");
        StsToolkit.runRunnable(new Runnable()
        {
            public void run()
            {
                rebuildModel(statusPanel);
                StsStatusPanel.disposeDialog();
            }
        });
    }

    /**
     * Set the string that represents the Z domain, generally time or depth.
     *
     * @params domain the byte representing Z (default is time)
     */
    private void setZDomainAndRebuild(byte domain)
    {
        setIsDepth(domain == TD_DEPTH);
        if(zDomain == domain) return;
        setZDomain(domain);
        if(!isModelBuilt()) return;
        rebuildModel();
    }

    public void setZDomain(byte domain)
    {
        setIsDepth(domain == TD_DEPTH);
        if(zDomain == domain) return;
        byte oldZDomain = zDomain;
        //        zDomainBean.setSelectedItem(StsParameters.TD_ALL_STRINGS[domain]);
        zDomain = domain;
        dbFieldChanged("zDomain", zDomain);

        StsWin3dFull[] parentWindows = currentModel.getParentWindows();
        if(parentWindows != null)      //for loading project that exited abnormally, prevents NPE
            for(StsWin3dFull parentWindow : parentWindows)
                parentWindow.setZDomain(domain);

        convertGridZ(oldZDomain);

        setBoundingBoxesZRanges(oldZDomain);
        resetCropVolume();
        if(!windowInitialized()) return;
        //This looks redundant (z slider adjusted above and resetCropVolume() called).  Tom 12/12/07
        rangeChanged();
        changeProjectView();
    }

    private void convertGridZ(byte prevZDomain)
    {
        displayBoundingBox.gridZ = convertZ(prevZDomain, displayBoundingBox.gridZ);
    }

    public boolean rebuildModel(StsStatusPanel statusPanel)
    {
        StsBuiltModel builtModel = (StsBuiltModel) currentModel.getCurrentObject(StsBuiltModel.class);
        builtModel.rebuildModel(currentModel, statusPanel);
        return true;
    }

    /**
     * Called when the project has changed between time and depth.
     * The new domain is defined by zDomain.
     */
    public void changeProjectView()
    {
        for(int i = 0; i < currentModel.viewPersistManager.families.length; i++)
        {
            StsWin3dBase[] windows = currentModel.getWindows(i);
            for(int w = 0; w < windows.length; w++)
            {
                StsWin3dBase window = windows[w];
                window.getCursor3d().resetDepthLabels(StsParameters.TD_ALL_STRINGS[zDomain]);
                StsView[] views = window.getDisplayedViews();
                for(int n = 0; n < views.length; n++)
                {
                    StsView view = views[n];
                    if(view instanceof StsView3d)
                    {
                        StsView3d view3d = (StsView3d) view;
                        float[] viewParams = view3d.getCenterAndViewParameters();
                        // these are the old values
                        float centerZ = viewParams[2];
                        float distance = viewParams[3];
                        float azimuth = viewParams[4];
                        float elevation = viewParams[5];
                        float zscale = viewParams[6];
                        float dz = distance * (float) Math.sin(elevation);
                        // compute the fractional height in the old domain and compute the height in the new domain
                        // then adjust the scaling from the old domain to the new domain
                        if(zDomain == TD_DEPTH)
                        {
                            float f = (centerZ - timeMin) / (timeMax - timeMin);
                            viewParams[2] = depthMin + f * (depthMax - depthMin);
                            zscale *= (timeMax - timeMin) / (depthMax - depthMin);
                            viewParams[6] = zscale;
                        }
                        else
                        {
                            float f = (centerZ - depthMin) / (depthMax - depthMin);
                            viewParams[2] = timeMin + f * (timeMax - timeMin);
                            zscale *= (depthMax - depthMin) / (timeMax - timeMin);
                            viewParams[6] = zscale;
                        }
                        view3d.setCenterAndViewParameters(viewParams);
                        view3d.glPanel3d.viewChanged();
                    }
                }
            }
        }
        //        currentModel.win3d.getCursor3d().rangeChanged();
        currentModel.resetAllSliderValues();
        currentModel.win3dDisplayAll();
    }

    public void setModelZDomainToCurrent()
    {
        setModelZDomain(zDomain);
    }

    public void setModelZDomain(byte zDomain)
    {
        if(modelZDomain == zDomain) return;
        modelZDomain = zDomain;
        dbFieldChanged("modelZDomain", modelZDomain);
    }

    private float convertZ(byte prevZDomain, float z)
    {
        if(prevZDomain == this.TD_TIME)
            return getDepthFromTimeFraction(z);
        else
            return getTimeFromDepthFraction(z);
    }

    public float getTimeFromDepth(float z)
    {
        if(timeMin == StsParameters.largeFloat) return 0.0f;
        float f = (z - depthMin) / (depthMax - depthMin);
        return timeMin + f * (timeMax - timeMin);
    }

    public float getDepthFromTime(float t)
    {
        if(timeMin == StsParameters.largeFloat) return 0.0f;
        float f = (t - timeMin) / (timeMax - timeMin);
        return depthMin + f * (depthMax - depthMin);
    }

    public float getTimeFromDepthFraction(float z)
    {
        if(timeMin == StsParameters.largeFloat) return 0.0f;
        float f = (z - depthMin) / (depthMax - depthMin);
        return timeMin + f * (timeMax - timeMin);
    }

    public float getDepthFromTimeFraction(float t)
    {
        if(timeMin == StsParameters.largeFloat) return 0.0f;
        float f = (t - timeMin) / (timeMax - timeMin);
        return depthMin + f * (depthMax - depthMin);
    }

    public boolean getHighResolution()
    {
        return highResolution;
    }

    public void setHighResolution(boolean value)
    {
        highResolution = value;
        //		setDisplayField("highResolution", highResolution);
        //		dbFieldChanged("highResolution", highResolution);
    }

    private float getZScaleMult(byte newDomain)
    {
        if(newDomain == TD_DEPTH)
        {
            return (depthMax - depthMin) / (timeMax - timeMin);
        }
        else
        {
            return (timeMax - timeMin) / (depthMax - depthMin);
        }
    }

    public boolean velocityUnitsChanged()
    {
        return velocityUnitsChanged(false);

    }

    public boolean velocityUnitsChanged(boolean silent)
    {
        if(!StsParameters.velocityUnitsChanged(defaultDistUnits, defaultTimeUnits, depthUnits, timeUnits))
            return false;
        if(!silent)
        {
            new StsMessage(currentModel.win3d, StsMessage.INFO,
                    "Velocity units being converted from " + StsParameters.getVelocityString(defaultDistUnits, defaultTimeUnits) + " to " +
                            StsParameters.getVelocityString(depthUnits, timeUnits));
        }
        return true;
    }

    public float convertVelocity(float velocity, StsModel model)
    {
        return StsParameters.convertVelocity(velocity, defaultDistUnits, defaultTimeUnits, depthUnits, timeUnits);
    }

    public float convertVelocityToProjectUnits(float velocity, String inputVelocityUnits)
    {
        return calculateVelScaleMultiplier(inputVelocityUnits) * velocity;
    }

    public boolean timeUnitsChanged()
    {
        if(!StsParameters.timeUnitsChanged(defaultTimeUnits, timeUnits)) return false;
        new StsMessage(currentModel.win3d, StsMessage.INFO,
                "Time units being converted from " + StsParameters.getTimeString(defaultTimeUnits) + " to " +
                        StsParameters.getTimeString(timeUnits));
        return true;
    }

    public float convertTime(float time, StsModel model)
    {
        return StsParameters.convertTime(time, defaultTimeUnits, timeUnits);
    }

    public boolean depthUnitsChanged()
    {
        if(!StsParameters.distanceUnitsChanged(defaultDistUnits, depthUnits)) return false;
        new StsMessage(currentModel.win3d, StsMessage.INFO,
                "Depth units being converted from " + StsParameters.getDepthString(defaultDistUnits) + " to " +
                        StsParameters.getDepthString(depthUnits));
        return true;
    }

    public float convertDepth(float depth, StsModel model)
    {
        return StsParameters.convertDistance(depth, defaultDistUnits, depthUnits);
    }

    private void resetWindowTimeOrDepth(StsWin3dBase window, byte newZDomain)
    {
    }

    private boolean windowInitialized()
    {
        return currentModel != null && currentModel.win3d != null && currentModel.win3d.getCursor3d() != null;
    }

    public byte getZDomainSupported()
    {
        return zDomainSupported;
    }

    public String[] getZDomainSupportedStrings()
    {
        return StsParameters.getSupportedDomainStrings(zDomainSupported);
    }

    private void setZDomainSupported(byte zDomainSupported)
    {
        if(this.zDomainSupported == zDomainSupported) return;
        this.zDomainSupported = zDomainSupported;
        zDomainSupportedBean.setValueObject(StsParameters.TD_ALL_STRINGS[zDomainSupported]);
        StsWin3dFull[] parentWindows = currentModel.getParentWindows();
        if(parentWindows != null)  //for loading project that exited abnormally, prevents NPE
            for(StsWin3dFull parentWindow : parentWindows)
                parentWindow.setZDomainSupported(zDomainSupported);
        dbFieldChanged("zDomainSupported", this.zDomainSupported);
    }

    /**
     * zDomainSupported indicates which domains can be handled: time, depth, or both.
     * zDomain indicates which domain is the current state: time or depth.
     * newDomain is the domain(s) which are supported by this object being added.
     * If not initialized (ZDOMAIN_NONE), then the supportedZDomain is newDomain
     * and zDomain is also newDomain.
     * If newDomain and zDomain agree, set supportedZDomain tothe same.
     * If newDomain is both and zDomain is one or the other,
     * then zDomain remains the same but supportedZDomain is set to zDomain.
     */
    static boolean retval = false;

    public boolean checkSetZDomain(byte newZDomainSupported)
    {
        final byte d = newZDomainSupported;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                retval = doCheckSetZDomain(d);
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
        return retval;
    }

    private boolean doCheckSetZDomain(byte newZDomainSupported)
    {
        if(newZDomainSupported == TD_TIME_DEPTH)
            return checkSetZDomain(newZDomainSupported, zDomain);
        else
            return checkSetZDomain(newZDomainSupported, newZDomainSupported);
    }

    /**
     * changes the supported zDomains for the project (TIME, DEPTH, or TIME_DEPTH), and
     * sets the current project zDomain (TIME or DEPTH).
     */

    public boolean checkSetZDomain(byte newZDomainSupported, byte preferZDomain)
    {
        if(debug)
        {
            System.out.println("Current: zDomain " + getZDomainString(zDomain) + " zDomainSupported " +
                    getZDomainString(zDomainSupported));
            System.out.println("   Input: zDomainSupported " + getZDomainString(newZDomainSupported) +
                    "   preferZDomain " +
                    getZDomainString(preferZDomain));
        }

        boolean ok = true;
        if(zDomainSupported == TD_NONE)
        {
            setZDomainSupported(newZDomainSupported);
        }
        else if(zDomainSupported != newZDomainSupported) // domains are different so we support all
        {
            setZDomainSupported(TD_TIME_DEPTH);
        }

        if(ok)
        {
            if(preferZDomain != TD_NONE)
            {
                setZDomainAndRebuild(preferZDomain);
            }
            else if(newZDomainSupported != TD_TIME_DEPTH)
            {
                setZDomainAndRebuild(newZDomainSupported);
            }
            else // should not occur
            {
                StsException.systemError("StsProject.checkZDomain() called with preferZDomain argument of " +
                        StsParameters.TD_ALL_STRINGS[preferZDomain] +
                        " must be either TIME or DEPTH.");
                setZDomainAndRebuild(TD_DEPTH);
            }
        }
        if(debug)
        {
            System.out.println("   New zDomain status: zDomain " + zDomain + " supportedZDomain " + zDomainSupported);
            System.out.println("   Status: " + ok);
        }
        //		zDomainBean.setValueObject(StsParameters.TD_ALL_STRINGS[zDomain]);
        //		zDomainBean.setEnabled(supportedZDomain == TD_TIME_DEPTH);
        //        resetProjectView();
        currentModel.refreshObjectPanel();
        return ok;
    }

    public String getZDomainString()
    {
        return StsParameters.TD_ALL_STRINGS[zDomain];
    }

    public String getZDomainString(byte zDomain)
    {
        return StsParameters.TD_ALL_STRINGS[zDomain];
    }

    public boolean isDepth()
    {
        return zDomain == TD_DEPTH;
    }

    public int getPointsZIndex()
    {
        if(zDomain == TD_DEPTH) // points are x,y,d,m and we want z index
        {
            return 2;
        }
        else // points are x,y,d,m,t and we want t index
        {
            return 4;
        }
    }

    public boolean isInProjectBounds(double x, double y)
    {
        float[] xy = getRotatedRelativeXYFromUnrotatedAbsoluteXY(x, y);
        return rotatedBoundingBox.isInsideXY(xy);
        /*
           double[] absMax = getAbsoluteXYCoordinates(getXMax(), getYMax());
           double[] absMin = getAbsoluteXYCoordinates(getXMin(), getYMin());
           if((x > absMax[0]) || (x < absMin[0]))
               return false;

           if((y > absMax[1]) || (y < absMin[1]))
               return false;

           return true;
       */
    }

    /*
          public String getLabelFormatAsString()
          {
              return labelFormat.toPattern();
          }

          public void setLabelFormat(DecimalFormat fmt)
          {
              labelFormat = fmt;
              return;
          }

          public void setLabelFormatAsString(String fmt)
          {
              labelFormat = new DecimalFormat(fmt);
              setDisplayField("labelFormatAsString", fmt);
              currentModel.win3dDisplayAll();
          }
      */
    public StsDisplayBoundingBox getDisplayBoundingBox()
    {
        return displayBoundingBox;
    }

    public StsRotatedGridBoundingBox getRotatedBoundingBox()
    {
        return rotatedBoundingBox;
    }

    public StsRotatedGridBoundingBox getZDomainRotatedBoundingBox(byte zDomain)
    {
        if(zDomain == TD_TIME && this.supportsTime())
            return getTimeRotatedBoundingBox();
        else if(zDomain == TD_DEPTH && supportsDepth())
            return getDepthRotatedBoundingBox();
        else
            return null;
    }

    public StsRotatedGridBoundingBox getTimeRotatedBoundingBox()
    {
        StsRotatedGridBoundingBox boundingBox = new StsRotatedGridBoundingBox(rotatedBoundingBox, false);
        boundingBox.resetZRange(timeMin, timeMax, timeInc);
        return boundingBox;
    }

    public StsRotatedGridBoundingBox getDepthRotatedBoundingBox()
    {
        StsRotatedGridBoundingBox boundingBox = new StsRotatedGridBoundingBox(rotatedBoundingBox, false);
        boundingBox.resetZRange(depthMin, depthMax, depthInc);
        return boundingBox;
    }

    public StsBoundingBox getUnrotatedBoundingBox()
    {
        return unrotatedBoundingBox;
    }

    public StsCropVolume getCropVolume()
    {
        return cropVolume;
    }

    public void resetCropVolume()
    {
        cropVolume.initialize();
        currentModel.instanceChange(cropVolume, "cropVolumeChanged");
    }

    /**
     * Is the value in between the supplied minimum and maximum. If not, display user error message.
     *
     * @params value the test value
     * @params minValue the minimum range value
     * @params maxValue the maximum range value
     * @params message the error message if out of range
     * @returns true if in range, false if outside range
     */
    private boolean valueOK(float value, float minValue, float maxValue, String message)
    {
        if(StsMath.betweenInclusive(value, minValue, maxValue))
        {
            return true;
        }
        new StsMessage(null, StsMessage.ERROR, message);
        return false;
    }

    /**
     * Get the location of the data directory as a Universal Resource Locator
     *
     * @returns URL string to data directory
     */
    public String getDataDirURLName()
    {
        if(Main.isWebStart)
        {
            java.net.URL codeBaseURL = JNLPUtilities.getBasicService().getCodeBase();
            if(codeBaseURL == null)
            {
                return new String("");
            }
            else
            {
                return codeBaseURL.toString();
            }
        }
        else
        {
            return new String("file:" + getRootDirString());
        }
    }

    /**
     * Make the Project root directory and the default data, binary and model folders under it.
     *
     * @params dir the parent directory for the Project as selected by the user
     */
    public void setRootDirectory(File dir)
    {
        try
        {
            if(dir == null)
            {
                rootDirectory = new File(".");
            }
            else if(!dir.isDirectory())
            {
                rootDirectory = new File(dir.getParent());
            }
            else
            {
                try
                {
                    rootDirectory = new File(dir.getCanonicalPath());
                }
                catch(IOException e)
                {
                    rootDirectory = new File(dir.getAbsolutePath());
                }
            }
            makeSubdir(rootDirectory, DATA_FOLDER);
            makeSubdir(rootDirectory, BINARY_FOLDER);
            makeSubdir(rootDirectory, MODEL_FOLDER);
            makeSubdir(rootDirectory, ARCHIVE_FOLDER);
            makeSubdir(rootDirectory, MEDIA_FOLDER);
            //			makeSubdir(Sts_3DMODELS_FOLDER);
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.setRoodDirectory() failed.",
                    e, StsException.WARNING);
        }
    }

    /**
     * Get the fully qualified path of the project root directory. If not set, it will be defaulted to user.dir
     *
     * @returns the pathname to the project root directory
     */
    public File getRootDirectory()
    {
        if(rootDirectory != null)
        {
            return rootDirectory;
        }
        else
        {
            return new File(System.getProperty("user.dir") + File.separator + "S2SCache" + File.separator);
        }
    }

    /**
     * Get the fully qualified path of the project model directory. If not set, it will be defaulted to user.dir
     *
     * @returns the pathname the project model directory
     */
    public File getModelDirectory()
    {
        if(modelDirectory != null)
        {
            return modelDirectory;
        }
        else
        {
            return new File(System.getProperty("user.dir") + File.separator + "S2SCache" + File.separator);
        }
    }

    /**
     * Get the fully qualified path of the data directory
     *
     * @returns the pathname to the data directory
     */
    public File getDataDirectory()
    {
        return dataDirectory;
    }

    /**
     * Get the fully qualified path of the binary directory
     *
     * @returns the pathname to the binary directory
     */
    public File getBinaryDirectory()
    {
        return binaryDirectory;
    }

    /**
     * Get the fully qualified path of the 3D Models directory (not used)
     *
     * @returns the pathname to the 3D model directory
     */
    public File getExportModelsDirectory()
    {
        return exportModelsDirectory;
    }

    /**
     * Get the fully qualified path of the project root directory as a String.
     *
     * @returns the pathname to the project root directory
     */
    public String getRootDirString()
    {
        return rootDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the data directory as a String.
     *
     * @returns the pathname to the data directory
     */
    public String getDataDirString()
    {
        return dataDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the model directory as a String
     *
     * @returns the pathname to the model directory
     */
    public String getModelDirString()
    {
        return modelDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the binary directory as a String
     *
     * @returns the pathname to the binary directory
     */
    public String getBinaryDirString()
    {
        return binaryDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the archive directory as a String
     *
     * @returns the pathname to the archvie directory
     */
    public String getArchiveDirString()
    {
        return archiveDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the media directory as a String
     *
     * @returns the pathname to the media directory
     */
    public String getMediaDirString()
    {
        return mediaDirectory.getPath() + File.separator;
    }

    /**
     * Get the fully qualified path of the 3D Models directory as a String (not used)
     *
     * @returns the pathname to the 3D model directory
     */
    public String getExportModelsDirString()
    {
        return exportModelsDirectory.getPath() + File.separator;
    }

    /**
     * Set the fully qualified path of the data directory. This is the location
     * where all cached data is stored
     *
     * @params dir set the data directory to this
     */
    public void setDataDirectory(File dir)
    {
        if(dir == null)
        {
            dataDirectory = new File("." + File.separator);
        }
        else
        {
            dataDirectory = dir;
        }
    }

    /**
     * Get the fully qualified path of the data directory. Stored at the
     * projects level, above binary and model directories.
     *
     * @returns the pathname to the data directory
     */
    public String getDataFullDirString()
    {
        String dataDirString = getDataDirString();
        if(dataDirString.equals("." + File.separator))
        {
            return getRootDirString();
        }
        return getRootDirString() + File.separator + dataDirString;
    }

    /**
     * Set the fully qualified path of the model directory. The model directory is where the
     * database files are stored.
     *
     * @params dir set the model directory to
     */
    public void setModelDirectory(File dir)
    {
        if(dir == null)
        {
            modelDirectory = new File("." + File.separator);
        }
        else
        {
            modelDirectory = dir;
        }
    }

    /**
     * Get the fully qualified path of the model directory
     *
     * @return pathname the model directory
     */
    public String getModelFullDirString()
    {
        String modelDirString = getModelDirString();
        if(modelDirString.equals("." + File.separator))
        {
            return getRootDirString();
        }
        return getRootDirString() + File.separator + modelDirString;
    }

    /**
     * Set the fully qualified path of the binary directory. The binary directory is where all
     * cached binary files are stored.
     *
     * @parm pathname set the binary directory to
     */
    public void setBinaryDirectory(File dir)
    {
        if(dir == null)
        {
            binaryDirectory = new File("." + File.separator);
        }
        else
        {
            binaryDirectory = dir;
        }
    }

    /**
     * Get the fully qualified path of the binary directory
     *
     * @return pathname for the binary directory
     */
    public String getBinaryFullDirString()
    {
        String binaryDirString = getBinaryDirString();

        if(binaryDirString.equals("." + File.separator))
            return getRootDirString();
        return getRootDirString() + File.separator + binaryDirString;
    }

    /**
     * Get the fully qualified path of the model directory (not used)
     *
     * @return pathname the model directory
     */
    public String getExportModelsFullDirString()
    {
        String exportModelsDirString = getExportModelsDirString();
        if(exportModelsDirString.equals("." + File.separator))
        {
            return getRootDirString();
        }
        return getRootDirString() + File.separator + exportModelsDirString;
    }

    /**
     * Set the name of the current Project
     *
     * @param name the new Project name
     */
    public void setName(String name)
    {
        this.name = name;
    }

    /**
     * Get the Project name
     *
     * @return name of the Project
     */
    public String getName()
    {
        return name;
    }

    /** adjust the unrotatedBoundingBox so it just bounds the rotatedBoundingBox */
    public void adjustUnrotatedGrid()
    {
        unrotatedBoundingBox.addRotatedBoundingBox(rotatedBoundingBox);
    }

    /**
     * Adjust boundingBox to fit the displayBoundingBox if angle has not been set yet;
     * otherwise, don't make changes to boundingBox.
     */
    public void adjustBoundingBoxToDisplayBoundingBox()
    {
        if(!rotatedBoundingBox.getAngleSet())
        {
            rotatedBoundingBox.initialize(displayBoundingBox);
            //		if(currentModel != null && currentModel.win3d != null) currentModel.resetCursorPanel();
        }
    }

    /**
     * Starting with current boundingBox, adjust x, y, and z ranges so they are "nice".
     *
     * @parameters model the current StsModel object
     */
    public void adjustDisplayBoundingBoxXYRange()
    {
        displayBoundingBox.adjustBoundingBoxXYRange(unrotatedBoundingBox);

        if(currentModel != null)
        {
            //            StsChangeCmd cmd = new StsChangeCmd(model, this, "project", false);
            //            currentModel.addTransactionCmd("adjustBoundingBox", cmd);

            // SAJ - Conflict with need to maintain cursor position on new data load
            //            if(currentModel.win3d != null) currentModel.resetCursorPanel();
        }
    }

    /*
       private void initializeBoundingBoxesZRange()
       {
        float zMin = unrotatedBoundingBox.zMin;
        float zMax = unrotatedBoundingBox.zMax;
        double[] zScale = StsMath.niceScale(zMin, zMax, approxNumZGridIncrements, true);
        zMin = (float) (zScale[0]);
        zMax = (float) (zScale[1]);
        zDomain = TD_DEPTH;
        setBoundingBoxesZRanges(zMin, zMax, 0.0f, true);
       }
       */
    private void setBoundingBoxesZRanges(byte prevZDomain)
    {
        if(!resetBoundingBoxesZRanges()) return;

        StsWin3dFull[] parentWindows = currentModel.getParentWindows();
        if(parentWindows == null) return;
        for(int n = 0; n < parentWindows.length; n++)
        {
            if(parentWindows[n] == null) continue;

            // Must re-initialize the axis for Cursor views from time to depth and visa-versa
            StsView[] views = parentWindows[n].getDisplayedViews();
            for(int i = 0; i < views.length; i++)
            {
                if(views[i] instanceof StsViewCursor)
                    ((StsViewCursor) views[i]).initializeAxisRanges();
            }

            // Reconfigure the cursor
            StsCursor3d cursor3d = parentWindows[n].getCursor3d();
            float z = cursor3d.getCurrentDirCoordinate(StsCursor3d.ZDIR);
            z = convertZ(prevZDomain, z);
            parentWindows[n].adjustSlider(StsCursor3d.ZDIR, z);
            cursor3d.zRangeChanged();
            currentModel.viewPersistManager.families[n].adjustCursor(StsCursor3d.ZDIR, z);
        }
    }

    private void setBoundingBoxesZRanges()
    {
        if(!resetBoundingBoxesZRanges()) return;

        StsWin3dFull[] parentWindows = currentModel.getParentWindows();
        if(parentWindows == null) return;
        for(int n = 0; n < parentWindows.length; n++)
        {
            if(parentWindows[n] == null) continue;
            StsCursor3d cursor3d = parentWindows[n].getCursor3d();
            float z = cursor3d.getCurrentDirCoordinate(StsCursor3d.ZDIR);
            //			z = convertZ(prevZDomain, z);
            parentWindows[n].adjustSlider(StsCursor3d.ZDIR, z);
            cursor3d.zRangeChanged();
            currentModel.viewPersistManager.families[n].adjustCursor(StsCursor3d.ZDIR, z);
        }
    }

    /**
     * apply integral (nice) values of zMin, zMax, and zInc to various bounding boxes.
     * If adjustProjectZValues is false, don't adjust rotatedBoundingBox for example when data is seismic
     * and we want an exact cube fitting the seismic.
     */
    private void adjustBoundingBoxesZRange(byte zDomain, boolean adjustProjectZValues)
    {
        double[] range;
        if(zDomain == TD_TIME)
        {
            range = niceScale(timeMin, timeMax);
            float timeMin = (float) range[0];
            float timeMax = (float) range[1];
            float timeInc = (float) range[2];
            setBoundingBoxesZRanges(timeMin, timeMax, timeInc, adjustProjectZValues);
        }
        else
        {
            range = niceScale(depthMin, depthMax);
            float depthMin = (float) range[0];
            float depthMax = (float) range[1];
            float depthInc = (float) range[2];
            setBoundingBoxesZRanges(depthMin, depthMax, depthInc, adjustProjectZValues);
        }
    }

    private void setBoundingBoxesZRanges(float zMin, float zMax, float zInc, boolean adjustProjectZValues)
    {
        unrotatedBoundingBox.setZRange(zMin, zMax);
        displayBoundingBox.setZRange(zMin, zMax);
        if(adjustProjectZValues)
        {
            rotatedBoundingBox.setZRange(zMin, zMax, zInc);
        }
    }

    private boolean resetBoundingBoxesZRanges()
    {
        float zMin, zMax, zInc;
        if(isDepth)
        {
            zMin = depthMin;
            zMax = depthMax;
            zInc = depthInc;
        }
        else
        {
            zMin = timeMin;
            zMax = timeMax;
            zInc = timeInc;
        }
        if(zMin == StsParameters.largeFloat) return false;
        if(unrotatedBoundingBox != null) unrotatedBoundingBox.setZRange(zMin, zMax);
        if(displayBoundingBox != null) displayBoundingBox.setZRange(zMin, zMax);
        if(rotatedBoundingBox != null) rotatedBoundingBox.resetZRange(zMin, zMax, zInc);
        displayBoundingBox.gridZ = zMax;
        return true;
    }

    public void checkCreateZInc()
    {
        if(rotatedBoundingBox.zInc != 0.0f)
        {
            return;
        }
        double[] zScale = niceScale(rotatedBoundingBox.zMin, rotatedBoundingBox.zMax);
        float zInc = (float) zScale[2];
        if(rotatedBoundingBox.zInc == zInc)
        {
            return;
        }
        if(rotatedBoundingBox.zInc == zInc) return;
        rotatedBoundingBox.setZInc(zInc);
        rotatedBoundingBox.dbFieldChanged("zInc", zInc);
        //       currentModel.addMethodCmd(this, "setDBZInc", new Object[] { new Float(zInc) } );
    }

    /*
       public void checkCreateZInc()
       {
        if(rotatedBoundingBox.zInc != 0.0f) return;
       double[] zScale = StsMath.niceScale(rotatedBoundingBox.zMin, rotatedBoundingBox.zMax, approxNumZGridIncrements, true);
        float zInc = (float)zScale[2];
        rotatedBoundingBox.setZInc(zInc);
        currentModel.addMethodCmd(this, "setDBZInc", new Object[] { new Float(zInc) } );
       }
       */
    /*
      public void setDBZInc(Float zInc)
      {
          rotatedBoundingBox.setZInc(zInc.floatValue());
      }
  */
    /** Starting with current boundingBox, adjust x, y, and z grid lines so they are "nice".
     * @parameters model the current StsModel object
     * @see setValuesFromBoundingBox
     * @see StsMath.niceScale2
     */
    /*
       public void adjustBoundingBoxGridLines(StsModel model)
       {
        displayBoundingBox.adjustBoundingBoxGridLines(boundingBox);

        if(model != null)
        {
         StsChangeCmd cmd = new StsChangeCmd(model, this, "project", false);
         currentModel.addTransactionCmd("adjustBoundingBox", cmd);
        }
       }
       */
    /** Starting with current boundingBox, adjust x, y, and z ranges. Will test each of three pairs
     *  of minimum and maximum and will adjust the bounding box limits if required.
     * @params range X, Y and Z minimum and maximum pairs. In this order: XMin, XMax, YMin, YMax, ZMin, ZMax
     */
    /*
       public void adjustRange(float[] range)
       {
        rotatedBoundingBox.adjustRange(range);
       }
       */

    /**
     * Build the full URL for the associated filename
     *
     * @params filename relative binary filename
     */
    public URL getBinaryFileURL(String filename)
    {
        String urlName = "null";

        try
        {
            urlName = "file:" + getBinaryFullDirString() + "/" + filename;
            return new URL(urlName);
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.getBinaryFileURL() failed.\n" +
                    "Couldn't find url: " + urlName, e, StsException.WARNING);
            return null;
        }
    }

    /**
     * Convert the XY coordinates in an StsPoint to absolute coordinates
     *
     * @params point the StsPoint to convert to absolute
     * @returns absolute X and Y coordinate
     */
    public double[] getAbsoluteXYCoordinates(StsPoint point)
    {
        return getAbsoluteXYCoordinates(point.v);
    }

    /**
     * Convert XY coordinates to absolute coordinates
     *
     * @params x the relative X position
     * @params y the relative Y position
     * @returns absolute X and Y coordinate
     */
    public double[] getAbsoluteXYCoordinates(float[] xy)
    {
        return getAbsoluteXYCoordinates(xy[0], xy[1]);
    }

    /**
     * Convert XY coordinates to absolute coordinates
     *
     * @params x the relative X position
     * @params y the relative Y position
     * @returns absolute X and Y coordinate
     */
    public double[] getAbsoluteXYCoordinates(float x, float y)
    {
        return rotatedBoundingBox.getAbsoluteXY(x, y);
    }

    /**
     * Convert absolute XY coordinates to local coordinates
     *
     * @params x the relative X position
     * @params y the relative Y position
     * @returns relative X and Y coordinate
     */
    public float[] getRelativeXY(double x, double y)
    {
        return rotatedBoundingBox.getRelativeXY(x, y);
    }

    /**
     * Given xy in rotated relative coordinate system,
     * return unrotated xy relative to project origin.
     */
    public float[] getUnrotatedRelativeXYFromRotatedXY(float x, float y)
    {
        return rotatedBoundingBox.getUnrotatedRelativeXYFromRotatedXY(x, y);
    }

    public float[] getUnrotatedRelativeXYFromAbsXY(double x, double y)
    {
        return rotatedBoundingBox.getUnrotatedRelativeXYFromAbsXY(x, y);
    }

    public float[] getRotatedRelativeXYFromUnrotatedAbsoluteXY(double x, double y)
    {
        return rotatedBoundingBox.getRotatedRelativeXYFromUnrotatedAbsoluteXY(x, y);
    }

    public float[] getRotatedRelativeXYFromUnrotatedRelativeXY(float x, float y)
    {
        return rotatedBoundingBox.getRotatedRelativeXYFromUnrotatedRelativeXY(x, y);
    }

	public void rotatePoint(double[] point)
	{
		float[] xy = rotatedBoundingBox.getRotatedRelativeXYFromUnrotatedRelativeXY((float)point[0], (float)point[1]);
		point[0] = xy[0];
		point[1] = xy[1];
	}
    
    public boolean getDisplayCompass() { return displayCompass; }
    public void setDisplayCompass(boolean display)
    {
		if(displayCompass == display) return;
        displayCompass = display;
        this.dbFieldChanged("displayCompass", displayCompass);
        currentModel.win3dDisplay();
    }
    /*
       public StsPoint getLocalCoordinates(double xOrigin, double yOrigin, StsPoint point)
       {
        float[] xy = boundingBox.getRelativeXY(xOrigin + point.v[0], yOrigin + point.v[1]);
        StsPoint projectPoint = new StsPoint(point);
        projectPoint.v[0] = xy[0];
        projectPoint.v[1] = xy[1];
        return projectPoint;
       }
       */

    private void initializeBoundingBoxes()
    {
        rotatedBoundingBox.initializedXY = false;
        rotatedBoundingBox.initializedZ = false;
        initializeBoundingBoxesZRange();
        adjustDisplayBoundingBoxXYRange();
        adjustBoundingBoxToDisplayBoundingBox();
    }

    public boolean isInitialized()
    {
        return unrotatedBoundingBox.originSet;
    }

    public void initializeBoundingBoxesZRange()
    {
        float zMin = unrotatedBoundingBox.zMin;
        float zMax = unrotatedBoundingBox.zMax;
        double[] zScale = niceScale(zMin, zMax);
        zMin = (float) (zScale[0]);
        zMax = (float) (zScale[1]);
        displayBoundingBox.gridZ = zMax;
        unrotatedBoundingBox.setZRange(zMin, zMax);
        displayBoundingBox.setZRange(zMin, zMax);
        rotatedBoundingBox.setZRange(zMin, zMax);
        //        rotatedBoundingBox.setZRange(zMin, zMax, 0.0f);
    }

    public void adjustBoundingBoxes(boolean saveToDB, boolean adjustProjectZValues)
    {
        adjustBoundingBoxesZRange(zDomain, adjustProjectZValues);
        adjustDisplayBoundingBoxXYRange();
        adjustBoundingBoxToDisplayBoundingBox();
        if(debug)
        {
            StsException.systemDebug(this, "adjustBoundingBoxes", "boundingBoxes:");
            StsToolkit.objectToPrint(System.err, "  RotatedBoundingBox", this.rotatedBoundingBox, StsBoundingBox.class);
            StsToolkit.objectToPrint(System.err, "  unRotatedBoundingBox", unrotatedBoundingBox, StsBoundingBox.class);
            StsToolkit.objectToPrint(System.err, "  displayBoundingBox", displayBoundingBox, StsBoundingBox.class);
        }

        if(saveToDB)
        {
            dbFieldChanged("timeMin", timeMin);
            dbFieldChanged("timeMax", timeMax);
            dbFieldChanged("timeInc", timeInc);
            dbFieldChanged("depthMin", depthMin);
            dbFieldChanged("depthMax", depthMax);
            dbFieldChanged("depthInc", depthInc);
            currentModel.instanceChange(this.unrotatedBoundingBox, "adjustUnrotatedBoundingBox");
            currentModel.instanceChange(this.rotatedBoundingBox, "adjustRotatedBoundingBox");
            currentModel.instanceChange(this.displayBoundingBox, "adjustDisplayBoundingBox");
        }
        objectPanelChanged();
    }

    /**
     * Once a number of objects have been added to the model by a wizard for example,
     * we adjust the boundingBoxes and the display position and view.
     * This is more efficient than doing adjustments after every object is added as this is
     * more work in the database.
     */

    public void runCompleteLoading()
    {
        StsToolkit.runWaitOnEventThread(new Runnable()
        {
            public void run()
            {
                completeLoading();
            }
        });
    }

    private void completeLoading()
    {
        adjustBoundingBoxes(true, false); // adjust boundingBoxZRange if 2nd arg is true: do more testing
        StsWin3dFull[] parentWindows = currentModel.getParentWindows();
        for(StsWin3dFull parentWindow : parentWindows)
        {
            parentWindow.cursor3dPanel.gridCheckBoxSetVisibleAndSelected();
            parentWindow.cursor3d.initialize();
            parentWindow.cursor3d.resetInitialCursorPositions();
        }
        Iterator<StsView> familyWindowViewIterator = currentModel.getFamilyWindowViewIterator();
        while (familyWindowViewIterator.hasNext())
        {
            StsView view = familyWindowViewIterator.next();
            view.setDefaultView();
        }
        objectPanelChanged();
    }

    /**
     * add an unrotatedClass (e.g. well) to the list IF there currently is no angle set for
     * a rotated coordinate system which occurs when a seismic volume or surface is loaded
     */
    public void checkAddUnrotatedClass(Class c)
    {
        if(rotatedBoundingBox == null)
        {
            return;
        }
        if(rotatedBoundingBox.getAngleSet())
        {
            return;
        }

        StsClass unrotatedClass = currentModel.getCreateStsClass(c);
        unrotatedClasses = (StsClass[]) StsMath.arrayAddElementNoRepeat(unrotatedClasses, unrotatedClass, StsClass.class);
    }

    public boolean anyDependencies()
    {
        return true;
    }

    public void actionPerformed(ActionEvent e)
    {
        System.out.println("Project file initialized.");
    }

    /**
     * Make a directory if needed
     *
     * @params subdir directory to create if needed
     */
    static public boolean makeSubdir(File directory, String subdir)
    {
        if(subdir == null)
        {
            return false;
        }
        if(subdir.equals("."))
        {
            return true;
        }
        File dir = new File(directory, subdir);
        if(dir.exists())
        {
            if(!dir.isDirectory())
            {
                return false;
            }
            return true;
        }
        return dir.mkdir();
    }

    /** Display method for Project */
    public void display(StsGLPanel3d glPanel3d)
    {
        displayBoundingBox.display(glPanel3d, currentModel, rotatedBoundingBox.angle);
    }

    /**
     * Get the maximum Project dimensions
     *
     * @returns maximum projection distance in X, Y or Z
     */
    public float getMaxProjectDimension()
    {
        return rotatedBoundingBox.getDimensions();
    }

    /**
     * Vertical index is in increments of dZ with 0 at Z=0.0
     *
     * @param z value for which we want index just above
     * @return index above z
     */
    public int getIndexAbove(float z)
    {
        float sliceF = rotatedBoundingBox.getSliceCoor(z);
        return StsMath.below(sliceF);
    }

    /**
     * Vertical index is in increments of dZ with 0 at Z=0.0
     *
     * @param z value for which we want index just below
     * @return iBelow index below z
     */
    public int getIndexBelow(float z)
    {
        float sliceF = rotatedBoundingBox.getSliceCoor(z);
        return StsMath.above(sliceF);
    }


    public int getIndexAbove(float z, boolean isDepth)
    {
        float sliceF;
        if(isDepth)
            sliceF = (z - depthMin) / depthInc;
        else
            sliceF = (z - timeMin) / timeInc;
        return StsMath.floor(sliceF);
    }


    public int getIndexBelow(float z, boolean isDepth)
    {
        float sliceF;
        if(isDepth)
            sliceF = (z - depthMin) / depthInc;
        else
            sliceF = (z - timeMin) / timeInc;
        return StsMath.ceiling(sliceF);
    }


    /**
     * Get Z value at slice index
     *
     * @params index the slice index
     */
    public float getZAtIndex(int index)
    {
        return rotatedBoundingBox.getZCoor(index);
    }

    public float getZAtIndex(int index, boolean isDepth)
    {
        if(isDepth)
        {
            float z = (depthMin + index * depthInc);
            return StsMath.minMax(z, depthMin, depthMax);
        }
        else
        {
            float z = (timeMin + index * timeInc);
            return StsMath.minMax(z, timeMin, timeMax);
        }
    }

    public StsGridDefinition getGridDefinition()
    {
        return gridDefinition;
    }

    public void setGridDefinition(StsGridDefinition gridDefinition)
    {
        this.gridDefinition = gridDefinition;
        dbFieldChanged("gridDefinition", gridDefinition);
    }

    // TODO: Ask user to give us the xInc, yInc, etc Rather than setting it here unconditionally
    // TODO: as we may subsequently load an object which has row col numbering
    public void checkCursor3d()
    {
        if(rotatedBoundingBox.rowNumMin == StsParameters.nullValue)
        {
            currentModel.getCursor3d().setIsGridCoordinates(false);
            double[] range = niceScale(rotatedBoundingBox.xMin, rotatedBoundingBox.xMax);
            rotatedBoundingBox.xInc = (float) range[2];
            range = niceScale(rotatedBoundingBox.yMin, rotatedBoundingBox.yMax);
            rotatedBoundingBox.yInc = (float) range[2];
        }
    }

    public StsSeismicVolume constructVelocityVolume(StsSeismicVolume inputVelocityVolume,
                                                    float topTimeDatum, float topDepthDatum, float minVelocity, float maxVelocity, double scaleMultiplier, boolean useSonic, boolean useVelf, String[] velfList, StsProgressPanel panel)
    {
        panel.appendLine("Velocity volume construction is unimplemented.");
        return null;
    }

    public StsSeismicVelocityModel constructVelocityModelSV(float topTimeDatum, float topDepthDatum, float minVelocity, float maxVelocity, double scaleMultiplier, double newTimeInc, StsProgressPanel panel)
    {
        if(velocityModel != null) velocityModel.delete();
        try
        {
            velocityModel = new StsSeismicVelocityModel(topTimeDatum, topDepthDatum, minVelocity, maxVelocity,
                    scaleMultiplier, newTimeInc, panel);
            if(velocityModel != null)
            {
                panel.appendLine("Velocity model construction is complete.");
                fieldChanged("velocityModel", velocityModel);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.constructVelocityModel() failed.", e, StsException.WARNING);
            return null;
        }

        return velocityModel;

    }

    public StsSeismicVelocityModel constructVelocityModel(StsModelSurface[] surfacesWithMarkers, StsSeismicVolume inputVelocityVolume, float tMin,
														  float zMin, float minVel, float maxVel, double scaleMultiplier, float newTimeInc,
														  float[] oneWayConstantIntervalVelocities, float oneWayConstantBottomIntervalVelocity, boolean useWellControl,

														  float markerFactor, int gridType, StsProgressPanel panel)
    {
        try
        {
            if(velocityModel != null) velocityModel.delete();
            velocityModel = new StsSeismicVelocityModel(surfacesWithMarkers, inputVelocityVolume, tMin, zMin, minVel, maxVel,
                    scaleMultiplier, newTimeInc, oneWayConstantIntervalVelocities, oneWayConstantBottomIntervalVelocity, useWellControl, markerFactor, gridType, panel);

            panel.appendLine("Velocity model construction is complete.");
            fieldChanged("velocityModel", velocityModel);
            //panel.appendLine("Adjusting project Z ranges");
            //adjustZRange(velocityModel);
            return velocityModel;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.constructVelocityModel() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /*
          private void adjustZRange(StsSeismicVelocityModel velocityModel)
          {
              StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
              float timeMin = velocityModel.timeMin;
              float timeMax = velocityModel.timeMax;
          }
      */
    public StsSeismicVelocityModel getSeismicVelocityModel()
    {
        return velocityModel;
    }
    /*
      public StsBasemap getBasemap() { return basemap; }

      public void setBasemapTextureObject(StsCursor3dDisplayable textureObject)
      {
          basemap.setBasemapTextureObject(textureObject);
      }
  */
    // Methods for StsTreeObjectI interface

    public Object[] getChildren()
    {
        return new Object[0];
    }

    /**
     * Get the object panel. Object panel is where persistent data object properties
     * are exposed to the user
     *
     * @return object panel
     */
    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null)
            objectPanel = StsObjectPanel.constructor(this, true);
        //        else
        //			objectPanel.setViewObject(this);
        return objectPanel;
    }

    public void deleteObjectPanel()
    {
        objectPanel = null;
    }

    /**
     * Is the Object Panel tree selected.
     *
     * @return true if selected
     */
    public void treeObjectSelected()
    {
    }

    /** Getters and Setters for all the default settings and user preferences. */
    /*
       public boolean getDefaultIsVisible()  { return defaultIsVisible; }
       public void setDefaultIsVisible(boolean value) { defaultIsVisible = value; fieldChanged("defaultIsVisible", defaultIsVisible); }
       public boolean getDefaultShowGrid()  { return defaultShowGrid; }
       public void setDefaultShowGrid(boolean value) { defaultShowGrid = value;  fieldChanged("defaultShowGrid", defaultShowGrid); }
       public boolean getDefaultShowLabels()  { return defaultShowLabels; }
       public void setDefaultShowLabels(boolean value) { defaultShowLabels = value;  fieldChanged("defaultShowLabels", defaultShowLabels); }
       public Color getDefaultBackgroundColor()  { return defaultBackgroundColor.getColor(); }
       public void setDefaultBackgroundColor(Color color) { defaultBackgroundColor.setBeachballColors(color);  fieldChanged("defaultBackgroundColor", defaultBackgroundColor); }
       public Color getDefaultGridColor()  { return defaultGridColor.getColor(); }
       public void setDefaultGridColor(Color color) { defaultGridColor.setBeachballColors(color);  fieldChanged("defaultGridColor", defaultGridColor); }
       public Color getDefaultTimingColor()  { return defaultTimingColor.getColor(); }
       public void setDefaultTimingColor(Color color) { defaultTimingColor.setBeachballColors(color);  fieldChanged("defaultTimingColor", defaultTimingColor); }
       public String getDefaultLabelFormatAsString() { return defaultLabelFormat; }
       public void setDefaultLabelFormatAsString(String value) { defaultLabelFormat = value;  fieldChanged("defaultLabelFormat", defaultLabelFormat); }
       */

    //JKF	public int index() { return -1; }
    //JKF	public void setIndex(int index) { }
    /*
      public StsDataCubeMemory getDataCubeMemory()
      {
          if(dataCubeMemory == null) dataCubeMemory = new StsDataCubeMemory();
          return dataCubeMemory;
      }
  */
    public StsBlocksMemoryManager getBlocksMemoryManager()
    {
        if(blocksMemoryManager == null)
        {
            blocksMemoryManager = new StsBlocksMemoryManager(1000); // jbw
            StsFileBlocks.blocksMemoryManager = blocksMemoryManager;
        }
        return blocksMemoryManager;
    }


    public void close()
    {
        if(blocksMemoryManager != null) blocksMemoryManager.clearAllBlocks();
        if(directoryAliasesFile != null) directoryAliasesFile.writeFile();
    }

    public String getDirectoryAlias(String directory)
    {
        if(directoryAliasesFile == null)
        {
            String fileDirectory = this.getRootDirString() + File.separator;
            String filename = "s2s.user.directoryAliases";
            directoryAliasesFile = new StsAsciiTokensFile(fileDirectory, filename, 2);
        }
        String[] matchToken = directoryAliasesFile.getMatchToToken(directory);
        if(matchToken == null) return null;
        return matchToken[1];
    }

    public void addDirectoryAlias(String directory, String alias)
    {
        if(directoryAliasesFile == null)
        {
            String fileDirectory = this.getRootDirString() + File.separator;
            String filename = "s2s.user.directoryAliases";
            directoryAliasesFile = new StsAsciiTokensFile(fileDirectory, filename, 2);
        }
        directoryAliasesFile.addLineTokens(new String[]{directory, alias});
    }

    public String getDirectoryPath(String directoryKeyName)
    {
        if(directoryTableFile == null)
        {
            String fileDirectory = this.getRootDirString();
            String filename = "s2s.user.directoryTable";
            directoryTableFile = new StsAsciiTokensFile(fileDirectory, filename, 2);
        }
        String[] matchToken = directoryTableFile.getMatchToToken(directoryKeyName);
        if(matchToken == null) return null;
        return matchToken[1];
    }

    public void setDirectoryPath(String keyName, String directoryPath)
    {
        if(keyName == null || directoryPath == null) return;
        if(directoryTableFile == null)
        {
            String fileDirectory = this.getRootDirString();
            String filename = "s2s.user.directoryTable";
            directoryTableFile = new StsAsciiTokensFile(fileDirectory, filename, 2);
        }
        directoryTableFile.setToken(keyName, directoryPath);
        directoryTableFile.writeFile();
    }

    static public void main(String[] args)
    {
        boolean ok;
        StsProject project = new StsProject();

        System.out.println("input seismic");
        project.zDomainSupported = StsProject.TD_NONE;
        project.setZDomainAndRebuild(StsProject.TD_NONE);
        ok = project.checkSetZDomain(StsProject.TD_TIME, StsProject.TD_TIME);

        System.out.println("input wells no td curves");
        project.zDomainSupported = StsProject.TD_NONE;
        project.setZDomainAndRebuild(StsProject.TD_NONE);
        ok = project.checkSetZDomain(StsProject.TD_DEPTH, StsProject.TD_DEPTH);

        System.out.println("input wells no td curves, no preferred domain");
        project.zDomainSupported = StsProject.TD_NONE;
        project.setZDomainAndRebuild(StsProject.TD_NONE);
        ok = project.checkSetZDomain(StsProject.TD_DEPTH, StsProject.TD_NONE);

        System.out.println("input wells with td curves");
        project.zDomainSupported = StsProject.TD_NONE;
        project.setZDomainAndRebuild(StsProject.TD_NONE);
        ok = project.checkSetZDomain(StsProject.TD_TIME_DEPTH, StsProject.TD_DEPTH);

        System.out.println("we have wells with no tds; input seismic - should fail");
        project.zDomainSupported = StsProject.TD_DEPTH;
        project.setZDomainAndRebuild(StsProject.TD_DEPTH);
        ok = project.checkSetZDomain(StsProject.TD_TIME);

        System.out.println("we have seismic (no td model); input wells without tds - should fail");
        project.zDomainSupported = StsProject.TD_TIME;
        project.setZDomainAndRebuild(StsProject.TD_TIME);
        ok = project.checkSetZDomain(StsProject.TD_DEPTH);

        System.out.println("we have seismic (no td model); input wells with tds - ok; leave in time");
        project.zDomainSupported = StsProject.TD_TIME;
        project.setZDomainAndRebuild(StsProject.TD_TIME);
        ok = project.checkSetZDomain(StsProject.TD_TIME_DEPTH);

        System.out.println(
                "we have seismic with td model; input wells with tds - ok; leave in time, supports time/depth");
        project.zDomainSupported = StsProject.TD_TIME_DEPTH;
        project.setZDomainAndRebuild(StsProject.TD_TIME);
        ok = project.checkSetZDomain(StsProject.TD_TIME_DEPTH);
    }

    public boolean getIsPerspective()
    {
        return isPerspective;
    }

    public void setIsPerspective(boolean perspective)
    {
        if(isPerspective == perspective) return;
        isPerspective = perspective;
        isPerspectiveBean.setValue(perspective);

        if(currentModel == null) return;
        if(currentModel.win3d == null) return;

        currentModel.win3d.toggleProjection();
    }

    public boolean getDisplayContours()
    {
        return displayContours;
    }

    public void setDisplayContours(boolean displayContours)
    {
        this.displayContours = displayContours;
    }

    public float getContourInterval()
    {
        return contourInterval;
    }

    public void setContourInterval(float contourInterval)
    {
        this.contourInterval = contourInterval;
    }

	public int getnRows()
	{
		return nRows;
	}

	public void setnRows(int nRows)
	{
		this.nRows = nRows;
	}

	public int getnCols()
	{
		return nCols;
	}

	public void setnCols(int nCols)
	{
		this.nCols = nCols;
	}

	public float getRowGridSize()
	{
		return rowGridSize;
	}

	public void setRowGridSize(float rowGridSize)
	{
		this.rowGridSize = rowGridSize;
	}

	public float getColGridSize()
	{
		return colGridSize;
	}

	public void setColGridSize(float colGridSize)
	{
		this.colGridSize = colGridSize;
	}
}
/*
class PersistentStrings implements java.io.Serializable
{
    ArrayList stringsArrayList;

    PersistentStrings()
    {
        stringsArrayList = new ArrayList(0);
    }

    ArrayList getStringsArrayList()
    {
        return stringsArrayList;
    }

    static PersistentStrings read(String filename)
    {
        FileInputStream fis;
        try
        {
            fis = new FileInputStream(filename);
        }
        catch (FileNotFoundException e)
        {
            return new PersistentStrings();
        }
        try
        {
            ObjectInputStream ois = new ObjectInputStream(fis);
            PersistentStrings persistentStrings = (PersistentStrings) ois.readObject();
            fis.close();
            return persistentStrings;
        }
        catch (Exception e)
        {
            return new PersistentStrings();
        }
    }

    void write(String filename)
    {
        FileOutputStream fos;
        try
        {
            fos = new FileOutputStream(filename);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(this);
        }
        catch (Exception e)
        {
            StsException.outputException("StsProject.PersistentStrings() failed to write file: " + filename,
                                         e, StsException.WARNING);
            return;
        }
    }
*/
