/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie Class instantiated by the movie wizard.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

public class StsMovie extends StsMainObject implements StsTreeObjectI
 {
     /** Static movie state variable START = 0 */
     transient public static final int START = 0;
     /** Static movie state variable REVERSE = 1 */
     transient public static final int REVERSE = 1;
     /** Static movie state variable PAUSE = 2 */
     transient public static final int PAUSE = 2;
     /** Static movie state variable STOP = 3 */
     transient public static final int STOP = 3;
     /** Static movie state variable PLAY = 4 */
     transient public static final int PLAY = 4;
     /** Static movie state variable END = 1 */
     transient public static final int END = 5;

     transient public static final int ELEVATION = 0;
     transient public static final int AZIMUTH = 1;
     transient public static int INCREMENT = 0;
     transient public static int DECREMENT = 1;

     // Line Animation Types
     /** Static line movie orientation variable INLINE = 1 */
     public static final int INLINE = 1;
     /** Static line movie orientation variable XLINE = 0 */
     public static final int XLINE = 0;
     /** Static slice movie orientation variable ZDIR = 2 */
     public static final int ZDIR = 2;

     // PostStack3d Animation Types
     /** Static volume movie orientation variable VOXEL = 1 */
     transient public static final int VOXEL = 1;
     /** Static volume movie orientation variable POINTSET = 0 */
     transient public static final int POINTSET = 0;

     /** Static movie animation type line = 0 */
     transient public static final byte LINE = 0;
     /** Static movie animation type volume = 1 */
     transient public static final byte VOLUME = 1;

     transient public static final byte SINGLE = 0;
     transient public static final byte CUMINC = 1;
     transient public static final byte CUMULATIVE = 2;
     transient public static final byte[] cumTypes = {SINGLE, CUMINC, CUMULATIVE};
     transient public static final String[] cumTypeStrings = {"Single Frame", "Cumulative Over Increment", "Cumulative"};

     /** Static movie speed variable FULL_SPEED = -1 */
     transient public static final int FULL_SPEED = 50;
     transient private StsSliderBean slider = null;
     private int direction = INLINE;
     private float[] range = new float[2];
     private float increment = 1.0f;
     private int delay = FULL_SPEED;
     private boolean loop = false;
     private boolean cycleVolumes = false;
     private byte cumType = SINGLE;
     private int cumulativeIncrement = 1;
     private boolean elevAnim = false;
     private boolean azimAnim = false;
     private int elevStart = 0;
     private int azimStart = 0;
     private int elevInc = 5;
     private int azimInc = 5;
     private byte atype = LINE;
     private byte volumeType = -1;
     private StsPointList pointSet = null;

     transient private int currentAzimuth = 0;
     transient private int currentElevation = 0;
     transient private float current = 0.0f;
     transient private boolean isRunning = false;
     // transient private StsToolbar toolbar = null;
     transient private StsModel model = null;

 //	transient Runnable runDisplay = null;

     static StsObjectPanel objectPanel = null;
     /**
      * Instance of this object
      */
 //	transient public static StsMovie movie = null;

     static public StsFieldBean[] propertyFields = null;

     /**
      * Default constructor
      */
     public StsMovie()
     {
         /*
         runDisplay = new Runnable()
         {
             public void run()
             {
                 try
                 {
                     display();
                 }
                 catch(Exception e)
                 {
                     System.out.println("display error ");
                     e.printStackTrace();
                 }
             }
         };
         */
     }

     /**
      * constructor
      */
     public StsMovie(StsModel model)
     {
         super(false);
         this.model = model;
         this.setLoop(getMovieClass().getDefaultLoop());
         this.setCycleVolumes(getMovieClass().getDefaultCycleVolumes());
         /*
         runDisplay = new Runnable()
         {
             public void run()
             {
                 try
                 {
                     display();
                 }
                 catch(Exception e)
                 {
                     System.out.println("display error ");
                     e.printStackTrace();
                 }
             }
         };
         */

     }

     public boolean initialize(StsModel model)
     {
         try
         {
             this.model = model;

             /*
             if (runDisplay == null)
                 runDisplay = new Runnable()
                 {
                   public void run()
                   {
                     try
                     {
                       display();
                     }
                     catch(Exception e)
                    {
                       System.out.println("display error ");
                       e.printStackTrace();
                    }
                   }
                 };
             */
             return true;
         }
         catch(Exception e)
         {
             StsException.outputException("StsMovie.classInitialize() failed.", e, StsException.WARNING);
             return false;
         }
     }

     public void enableAnimation(boolean bool)
     {
         if(getType() == VOLUME)
             getPointSet().setIsAnimated(bool);
     }

     /**
      * Set the orientation that the movie is playing
      * @param dir INLINE=1, XLINE=0 or ZDIR=2
      */
     public void setDirection(int dir)
     {
         direction = dir;
     }

     public void setDirectionString(String dir)
     {
         direction = 1;
         if(dir.equalsIgnoreCase("XLINE") || dir.equalsIgnoreCase("Crossline"))
             direction = 0;
         else if(dir.equalsIgnoreCase("ZDIR") || dir.equalsIgnoreCase("Z") || dir.equalsIgnoreCase("Depth"))
             direction = 2;
     }

     /**
      * Set the range of frames to play. Range is specified in inline, xline or z slice ranges
      * @param dataRange min and max play range
      */
     public void setRange(float[] dataRange)
     {
         range[0] = dataRange[0];
         range[1] = dataRange[1];
     }

     public void setRange(float min, float max)
     {
         range[0] = min;
         range[1] = max;
     }

     public void setStartRange(float start)
     {
         float[] range = model.getProject().getRotatedBoundingBox().getBoundingBoxRangeData(direction);
         start = StsMath.minMax(start, range[0], range[1] - range[2]); // must be between min and max-inc
         setRange(start, range[1]);
     }

     public void setStopRange(float end)
     {
         float[] range = model.getProject().getRotatedBoundingBox().getBoundingBoxRangeData(direction);
         end = StsMath.minMax(end, range[0] + range[2], range[1]); // must be between min+inc and max
         setRange(range[0], end);
     }

     /**
      * Set the movie increment
      * @param inc frame increment (defaults to native step size)
      */
     public void setIncrement(float inc)
     {
         if(getType() == LINE)
         {
             float[] range = model.getProject().getRotatedBoundingBox().getBoundingBoxRangeData(direction);
             inc = Math.max(inc, range[2]);
         }
         increment = inc;
         dbFieldChanged("increment", increment);
     }

     /**
      * Set the delay between frames in milliseconds
      * @param delayMs delay(milliseconds)
      */
     public void setDelay(int delayMs)
     {
         delay = delayMs;
         dbFieldChanged("delay", delay);
     }

     /**
      * Set the starting elevation
      * @param angle(degrees)
      */
     public void setElevationStart(int angle)
     {
         elevStart = angle;
         dbFieldChanged("elevStart", elevStart);
     }

     /**
      * Set the starting azimuth
      * @param angle(degrees)
      */
     public void setAzimuthStart(int angle)
     {
         azimStart = angle;
         dbFieldChanged("azimStart", azimStart);
     }

     /**
      * Set the elevation increment
      * @param increment(degrees)
      */
     public void setElevationIncrement(int increment)
     {
         elevInc = increment;
         dbFieldChanged("elevInc", elevInc);
     }

     /**
      * Set the azimuth increment
      * @param increment(degrees)
      */
     public void setAzimuthIncrement(int increment)
     {
         azimInc = increment;
         dbFieldChanged("azimInc", azimInc);
     }

     /**
      * Set whether the movie includes elevation animation
      * @param animate true to elevation animate
      */
     public void setElevationAnimation(boolean animate)
     {
         elevAnim = animate;
         dbFieldChanged("elevAnim", elevAnim);
     }

     /**
      * Set whether the movie includes azimuth animation
      * @param animate true to animate azimuth
      */
     public void setAzimuthAnimation(boolean animate)
     {
         azimAnim = animate;
         dbFieldChanged("azimAnim", azimAnim);
     }

     /**
      * Set whether the movie loops at completion
      * @param lp true to loop
      */
     public void setLoop(boolean lp)
     {
         loop = lp;
         dbFieldChanged("loop", loop);
     }

     /**
      * Set whether the movie includes cycling volumes on each loop
      * @param cv true to cycle
      */
     public void setCycleVolumes(boolean cv)
     {
         cycleVolumes = cv;
         dbFieldChanged("cycleVolumes", cycleVolumes);
     }

     /**
      * Set whether the movie is cummlative or instance
      * @param cumm true to accumulate points
      */
     public void setCumulativeIncrement(int cumm)
     {
         cumulativeIncrement = cumm;
         dbFieldChanged("cumulativeIncrement", cumulativeIncrement);
     }

     /**
      * Set whether the movie is cummlative or instance
      */
     public void setCumTypeString(String type)
     {
         for(int i = 0; i < cumTypes.length; i++)
             if(type.equals(cumTypeStrings[i]))
             {
                 cumType = cumTypes[i];
                 dbFieldChanged("cumType", cumType);
                 return;
             }
     }

     public void setCumType(byte type)
     {
         if(type >= 0 && type < cumTypes.length)
            cumType = type;
     }


         /**
          * Turn on a set of frames based on user selection
          * Should only be temporary, until next setCurrentFrame or range
          */
       public void setFrameRange(int start, int end)
       {
           byte type = getType();
           int oldStart = getPointSet().getStartVolume();
           int oldEnd = getPointSet().getEndVolume();

           int cf = end;
           setType(CUMINC);
           start = start * (int)increment - 1;
           end = end * (int)increment - 1;
           getPointSet().setStartVolume(start);
           getPointSet().setEndVolume(end);

           StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { display(); }});
           /*
           if(javax.swing.SwingUtilities.isEventDispatchThread())
               display();
           else
           {
               try
               {
                   javax.swing.SwingUtilities.invokeAndWait(runDisplay);
               }
               catch(Exception e)
               {
                   StsException.systemError("StsMovie.setCurrentFrame() invokeAndWait() failed.");
               }
           }
           */
           setType(type);
       }

         /**
      * Set the current frame
      * @param cf any valid frame (aka slice number)
      */
     public void setCurrentFrame(int cf)
     {
 //        System.out.println("Set frame= " + cf);
         if(getType() == StsMovie.VOLUME)
         {
             if(cumType == CUMINC)
             {
                 if(cf > 0)
                     cf = cf * (int)increment - 1;
                 int start = cf - cumulativeIncrement + 1;
                 if(start < 0)
                     start = 0;
                 getPointSet().setStartVolume(start);
                 getPointSet().setEndVolume(cf);
             }
             else if(cumType == CUMULATIVE)
             {
                 if(cf > 0)
                     cf = cf * (int)increment - 1;
                 getPointSet().setStartVolume(0);
                 getPointSet().setNextCummlativeVolume((int)cf);
             }
             else
             {
                 if(cf > 0)
                     cf = cf * (int)increment - 1;
                 getPointSet().setCurrentVolume((int)cf);
             }
         }
         else
         {
             int volNumber = cf / getNumberSlices();
             int sliceNumber = cf % getNumberSlices();
             float zVal = getRange()[0] + (sliceNumber * getIncrement());
 //            System.out.println("PostStack3d number:" + volNumber + " Slice number:" + sliceNumber + " Z Value:" + zVal);
             if(cycleVolumes)
                 setVolume(volNumber);

             slider.setValue(zVal);
             current = sliceNumber;
         }
         StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { display(); }});
         /*
         if(javax.swing.SwingUtilities.isEventDispatchThread())
             display();
         else
         {
             try
             {
                 javax.swing.SwingUtilities.invokeAndWait(runDisplay);
             }
             catch(Exception e)
             {
                 StsException.systemError("StsMovie.setCurrentFrame() invokeAndWait() failed.");
             }
         }
         */
     }

     public void display()
     {
         model.win3dDisplayAll();
     }

     public void setFrame(int cf)
     {
         if(getType() == StsMovie.VOLUME)
         {
             if(cumType == CUMINC)
             {
                 int start = cf - cumulativeIncrement + 1;
                 if(start < 0)
                     start = 0;
                 getPointSet().setStartVolume(start);
                 getPointSet().setEndVolume(cf);
             }
             else if(cumType == CUMULATIVE)
                 getPointSet().setNextCummlativeVolume((int)cf);
             else
                 getPointSet().setCurrentVolume((int)cf);
         }
         current = cf;
         StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { display(); }});
         /*
         if(javax.swing.SwingUtilities.isEventDispatchThread())
             display();
         else
         {
             try
             {
                 javax.swing.SwingUtilities.invokeAndWait(runDisplay);
             }
             catch(Exception e)
             {
                 StsException.systemError("StsMovie.setCurrentFrame() invokeAndWait() failed.");
             }
         }
         */
     }

     public void initializeEnvironment()
     {
         prepareCursor();
         slider = prepareSlider();
     }

     /**
      * Set the current azimuth
      * @return azimuth
      */
     public void setCurrentAzimuth(int azimuth)
     {
         currentAzimuth = azimuth;
     }

     /**
      * Set the current elevation
      * @return elevation
      */
     public void setCurrentElevation(int elevation)
     {
         currentElevation = elevation;
     }

     /**
      * Set whether the movie is currently being played or not
      * @param sp true if playing
      */
     public void setPlaying(boolean sp)
     {
         isRunning = sp;
     }

     public void setType(byte animationType)
     {
         atype = animationType;
     }

     public byte getType()
     {
         return atype;
     }

     public void setVolumeType(byte volumeAnimType)
     {
         volumeType = volumeAnimType;
     }

     public byte getVolumeType()
     {
         return volumeType;
     }

     public void setPointSet(StsPointList pointset)
     {
         pointSet = pointset;
     }

     public StsPointList getPointSet()
     {
         return pointSet;
     }

     /**
      * Set the toolbar associated with this movie
      * @param tb movie toolbar instance
      */
     /*
     public void setToolbar(StsToolbar tb)
     {
         toolbar = tb;
     }
     */
     /**
      * Get the direction of the movie
      * @return INLINE=1, XLINE=0 or ZDIR=2
      */
     public int getDirection()
     {
         return direction;
     }

     public String getDirectionString()
     {
         switch(direction)
         {
             case INLINE:
                 return "Y or InLine";
             case XLINE:
                 return "X or Crossline";
             case ZDIR:
                 return "T or Z";
             default:
                 return null;
         }
     }

     /**
      * Get the frame range of the movie. Range is specified in inline, xline or z slice ranges
      * @return min and max
      */
     public float[] getRange()
     {
         return range;
     }

     public float getStartRange()
     {
         return range[0];
     }

     public float getStopRange()
     {
         return range[1];
     }

     public int getNumberSlices()
     {
         return(int)((range[1] - range[0]) / increment);
     }

     /**
      * Get the frame increment. Defaulted to the slice increment.
      * @return frame increment
      */
     public float getIncrement()
     {
         return increment;
     }

     /**
      * Get the frame delay
      * @return delay in milliseconds
      */
     public int getDelay()
     {
         return delay;
     }

     /**
      * Get the loop. Specifies whether the movie should loop upon completion
      * @return true if set to loop
      */
     public boolean getLoop()
     {
         return loop;
     }

     /**
      * Get the volume cycle. Specifies whether the movie should cycle volumes upon completion
      * @return true if set to cycle
      */
     public boolean getCycleVolumes()
     {
         return cycleVolumes;
     }

     /**
      * Get the cummulative flag. Specifies whether the movie hould accumulate points or replace each frame
      * @return true if set to accumulate
      */
     public byte getCumulativeType()
     {
         return cumType;
     }

     /**
      * Set whether the movie is cummlative or instance as String
      */
     public String getCumTypeString()
     {
         return cumTypeStrings[cumType];
     }

     /**
      * Get the increment cummulative flag. Specifies whether the movie should accumulate points within each increment
      * @return true if set to accumulate
      */
     public int getCumulativeIncrement()
     {
         return cumulativeIncrement;
     }

     /**
      * Determine whether the movie includes elevation animation
      * @return true if set to animate
      */
     public boolean getElevationAnimation()
     {
         return elevAnim;
     }

     /**
      * Determine whether the movie includes azimuth animation
      * @return true if set to animate
      */
     public boolean getAzimuthAnimation()
     {
         return azimAnim;
     }

     /**
      * Get the elevation animation start angle
      * @return angle in degrees
      */
     public int getElevationStart()
     {
         return elevStart;
     }

     /**
      * Get the elevation animation increment angle
      * @return increment in degrees
      */
     public int getElevationIncrement()
     {
         return elevInc;
     }

     /**
      * Get the azimuith animation start angle
      * @return angle in degrees
      */
     public int getAzimuthStart()
     {
         return azimStart;
     }

     /**
      * Get the azimuth animation increment angle
      * @return increment in degrees
      */
     public int getAzimuthIncrement()
     {
         return azimInc;
     }

     /**
      * Get the current frame
      * @return frame number (aka slice number)
      */
     public float getCurrentFrame()
     {
         return current;
     }

     /**
      * Get the current azimuth
      * @return azimuth
      */
     public int getCurrentAzimuth()
     {
         return currentAzimuth;
     }

     /**
      * Get the current elevation
      * @return elevation
      */
     public int getCurrentElevation()
     {
         return currentElevation;
     }

     /**
      * Is the movie currently playing
      * @return true if movie is actively playing
      */
     public boolean isPlaying()
     {
         return isRunning;
     }

     /**
      * Get the toolbar instance associated with this movie.
      * @return toolbar instance
      */
     /*
     public StsToolbar getToolbar()
     {
        return toolbar;
     }
     */
     public int getNumberFrames()
     {
         if(getType() == VOLUME)
         {
             if(getVolumeType() == POINTSET)
             {
                 return (int)((float)getPointSet().getNumberVolumes()/getIncrement());
             }
             else
             {
                 // Not implemeted yet.
             }
         }
         else
         {
             int nVols = 1;
             if(getCycleVolumes())
                 nVols = ((StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class)).length;

             return nVols * getNumberSlices();
         }
         return 1;
     }

     /**
      * Movie selected on the Object tree
      */
     public void treeObjectSelected()
     {
         getMovieClass().selected(this);
         if(!currentModel.win3d.hasView(StsView3d.class))
         {
             new StsMessage(currentModel.win3d, StsMessage.ERROR, "Add a 3d view for movies.");
             return;
         }
         // currentModel.getGlPanel3d().checkAddView(StsView3d.class);
         currentModel.win3dDisplayAll();
     }

     static public StsMovieClass getMovieClass()
     {
         return(StsMovieClass)currentModel.getCreateStsClass(StsMovie.class);
     }

     public boolean anyDependencies()
     {
         return false;
     }

     public StsFieldBean[] getDisplayFields()
     {
         return null;
     }

     public StsFieldBean[] getPropertyFields()
     {
         if(propertyFields == null)
         {
             propertyFields = new StsFieldBean[]
             {
                 new StsStringFieldBean(StsMovie.class, "name", true, "Name:"),
                 new StsStringFieldBean(StsMovie.class, "directionString", false, "Direction:"),
                 new StsFloatFieldBean(StsMovie.class, "startRange", true, "Start At:"),
                 new StsFloatFieldBean(StsMovie.class, "stopRange", true, "End At:"),
                 new StsFloatFieldBean(StsMovie.class, "increment", true, "Increment:"),
                 new StsIntFieldBean(StsMovie.class, "delay", true, "Delay:"),
                 new StsBooleanFieldBean(StsMovie.class, "elevationAnimation", "Animated Elevation"),
                 new StsIntFieldBean(StsMovie.class, "elevationStart", true, "Start Elevation:"),
                 new StsIntFieldBean(StsMovie.class, "elevationIncrement", true, "Elevation Increment:"),
                 new StsBooleanFieldBean(StsMovie.class, "azimuthAnimation", "Animated Azimuth"),
                 new StsIntFieldBean(StsMovie.class, "azimuthStart", true, "Start Azimuth:"),
                 new StsIntFieldBean(StsMovie.class, "azimuthIncrement", true, "Azimuth Increment:"),
                 new StsBooleanFieldBean(StsMovie.class, "loop", "Loop"),
                 new StsBooleanFieldBean(StsMovie.class, "cycleVolumes", "Cycle Volumes"),
                 new StsStringFieldBean(StsMovie.class, "cumTypeString", "Type:"),
                 new StsIntFieldBean(StsMovie.class, "cumulativeIncrement", true, "Cumulative Increment:")
             };
         }
         return propertyFields;
     }

     public Object[] getChildren()
     {
         return new Object[0];
     }

     public StsObjectPanel getObjectPanel()
     {
         if(objectPanel == null)objectPanel = StsObjectPanel.constructor(this, true);
         return objectPanel;
     }

     /**
      * Change the volume to the next one in sequence.
      */
     public void nextVolume()
     {
         if(getVolumeType() == StsMovie.POINTSET) return;
         StsWin3dBase win3d = model.win3d;
         if(!win3d.hasView(StsView3d.class) && !win3d.hasView(StsViewCursor.class)) return;
         StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
         StsSeismicVolumeClass svc = (StsSeismicVolumeClass)model.getCreateStsClass(StsSeismicVolume.class);
         StsSeismicVolume currentSeismicVolume = svc.getCurrentSeismicVolume();
         for(int i = 0; i < volumes.length; i++)
         {
             if(currentSeismicVolume == volumes[i])
             {
                 if(i < (volumes.length - 1))
                     currentSeismicVolume = volumes[i + 1];
                 else
                     currentSeismicVolume = volumes[0];

                 svc.setCurrentObject(currentSeismicVolume);
                 StsComboBoxToolbar ctb = (StsComboBoxToolbar)model.win3d.getToolbarNamed(StsComboBoxToolbar.NAME);
                 if(ctb != null)
                     ctb.comboBoxSetItem(currentSeismicVolume);
                 return;
             }
         }
     }

     /**
      * Change the volume to the next one in sequence.
      */
     public void setVolume(int idx)
     {
         StsSeismicVolume[] volumes = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
         model.setCurrentObject(volumes[idx]);
     }

     /**
      * Initialize the particular movie slider and return it.
      * @return the slider corresponding to the movie definition.
      */
     public StsSliderBean prepareSlider()
     {
         StsSliderBean slider = null;
         if(getType() == StsMovie.VOLUME)
         {
             slider = null;
         }
         else
         {
             if(getDirection() == INLINE)
             {
                 slider = model.win3d.cursor3dPanel.getSliderY();
             }
             if(getDirection() == XLINE)
             {
                 slider = model.win3d.cursor3dPanel.getSliderX();
             }
             if(getDirection() == ZDIR)
             {
                 slider = model.win3d.cursor3dPanel.getSliderZ();
             }
             if(slider != null)
             {
                 if(!slider.isSelected())
                     slider.setSelected(true);
                 slider.setIncrement(getIncrement());
             }
         }
         return slider;
     }

     /**
      * Initialize the graphics cursor to the proper movie frame
      */
     public void prepareCursor()
     {
         if(getVolumeType() != StsMovie.POINTSET)
         {
             if(!model.win3d.isCursor3dDisplayed)
                 model.win3d.toggle3dCursor();
             if(model.win3d.cursor3dPanel != null)
                 model.win3d.cursor3dPanel.setSliderValues();
             if(!model.win3d.cursor3dPanel.gridChk.isSelected())
                 model.win3d.cursor3dPanel.gridChk.setSelected(true);
         }
     }

     /**
      * Increment or decrement the elevation angle.
      */
     public void setElevationAndAzimuth(int dir, int elevAzim)
     {
         int angle = 0, increment = 0, current = 0;

         StsView view = model.win3d.getCurrentView();
         StsView3d view3d = null;

         if(view instanceof StsView3d)
             view3d = (StsView3d)view;

         // Determine the current angle and increment
         switch(elevAzim)
         {
             case ELEVATION:
                 current = getCurrentElevation();
                 increment = getElevationIncrement();
                 break;
             case AZIMUTH:
                 current = getCurrentAzimuth();
                 increment = getAzimuthIncrement();
                 break;
             default:
                 return;
         }
         // Compute the angle
         if(dir == INCREMENT)
         {
             angle = current + increment;
             if(angle > 360)
                 angle -= 360;
         }
         else
         {
             angle = current - increment;
             if(angle < 0)
                 angle += 360;
         }
         // Change the model view
         switch(elevAzim)
         {
             case ELEVATION:
                 setCurrentElevation(angle);
                 if(view != null)view3d.setViewElevation(angle);
                 break;
             case AZIMUTH:
                 setCurrentAzimuth(angle);
                 if(view != null)view3d.setViewAzimuth(angle);
                 break;
             default:
                 return;
         }
     }

 }
