
/**
 * <p>Title: S2S Development</p>
 * <p>Description: Movie Class instantiated by the movie wizard.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */

package com.Sts.Types;

import com.Sts.UI.Toolbars.*;

public class StsMovie
{
    /** Static movie state variable START = 0 */
    public static final int START = 0;
    /** Static movie state variable REVERSE = 1 */
    public static final int REVERSE = 1;
    /** Static movie state variable PAUSE = 2 */
    public static final int PAUSE = 2;
    /** Static movie state variable STOP = 3 */
    public static final int STOP = 3;
    /** Static movie state variable PLAY = 4 */
    public static final int PLAY =  4;
    /** Static movie state variable END = 1 */
    public static final int END = 5;
    /** Static movie orientation variable INLINE = 1 */
    public static final int INLINE = 1;
    /** Static movie orientation variable XLINE = 0 */
    public static final int XLINE = 0;
    /** Static movie orientation variable ZDIR = 2 */
    public static final int ZDIR = 2;
    /** Static movie speed variable FULL_SPEED = -1 */
    public static final int FULL_SPEED = 50;

    private int direction = INLINE;
    private float[] range = new float[2];
    private float increment = 1.0f;
    private int delay = FULL_SPEED;
    private boolean loop = false;
    private boolean elevAnim = false;
    private boolean azimAnim = false;
    private int elevStart = 0;
    private int azimStart = 0;
    private int elevInc = 5;
    private int azimInc = 5;
    private int currentAzimuth = 0;
    private int currentElevation = 0;
    private float current = 0.0f;
    private boolean isRunning = false;
    private StsToolbar toolbar = null;
    /**
     * Instance of this object
     */
    public static StsMovie movie = null;
    /**
     * Default constructTraceAnalyzer
     */
    public void StsMovie()
    {

    }
    /**
     * Get the current movie object
     * @return movie object
     */
    static public StsMovie getInstance()
    {
        if(movie != null) return movie;
        movie = new StsMovie();
        return movie;
    }
    /**
     * Set the orientation that the movie is playing
     * @param dir INLINE=1, XLINE=0 or ZDIR=2
     */
    public void setDirection(int dir) { direction = dir; }
    /**
     * Set the range of frames to play. Range is specified in inline, xline or z slice ranges
     * @param dataRange min and max play range
     */
    public void setRange(float[] dataRange)
    {
        range[0] = dataRange[0];
        range[1] = dataRange[1];
    }
    /**
     * Set the movie increment
     * @param inc frame increment (defaults to native step size)
     */
    public void setIncrement(float inc) { increment = inc; }
    /**
     * Set the delay between frames in milliseconds
     * @param delayMs delay(milliseconds)
     */
    public void setDelay(int delayMs) { delay = delayMs; }
    /**
     * Set the starting elevation
     * @param angle(degrees)
     */
    public void setElevationStart(int angle) { elevStart = angle; }
    /**
     * Set the starting azimuth
     * @param angle(degrees)
     */
    public void setAzimuthStart(int angle) { azimStart = angle; }
    /**
     * Set the elevation increment
     * @param increment(degrees)
     */
    public void setElevationIncrement(int increment) { elevInc = increment; }
    /**
     * Set the azimuth increment
     * @param increment(degrees)
     */
    public void setAzimuthIncrement(int increment) { azimInc = increment; }
    /**
     * Set whether the movie includes elevation animation
     * @param animate true to elevation animate
     */
    public void setElevationAnimation(boolean animate) { elevAnim = animate; }
    /**
    * Set whether the movie includes azimuth animation
    * @param animate true to animate azimuth
    */
    public void setAzimuthAnimation(boolean animate) { azimAnim = animate; }
    /**
     * Set whether the movie loops at completion
     * @param lp true to loop
     */
    public void setLoop(boolean lp) { loop = lp; }
    /**
     * Set the current frame
     * @param cf any valid frame (aka slice number)
     */
    public void setCurrentFrame(float cf) { current = cf; }
    /**
     * Set the current azimuth
     * @return azimuth
     */
    public void setCurrentAzimuth(int azimuth) { currentAzimuth = azimuth; }
    /**
     * Set the current elevation
     * @return elevation
     */
    public void setCurrentElevation(int elevation) { currentElevation = elevation; }
    /**
     * Set whether the movie is currently being played or not
     * @param sp true if playing
     */
    public void setPlaying(boolean sp) { isRunning = sp; }
    /**
     * Set the toolbar associated with this movie
     * @param tb movie toolbar instance
     * @see StsMovieActionToolbar(StsModel, boolean)
     */
    public void setToolbar(StsToolbar tb) { toolbar = tb; }

    /**
     * Get the direction of the movie
     * @return INLINE=1, XLINE=0 or ZDIR=2
     */
    public int getDirection() { return direction; }
    /**
     * Get the frame range of the movie. Range is specified in inline, xline or z slice ranges
     * @return min and max
     */
    public float[] getRange() { return range; }
    /**
     * Get the frame increment. Defaulted to the slice increment.
     * @return frame increment
     */
    public float getIncrement() { return increment; }
    /**
     * Get the frame delay
     * @return delay in milliseconds
     */
    public int getDelay() { return delay; }
    /**
     * Get the loop. Specifies whether the movie should loop upon completion
     * @return true if set to loop
     */
    public boolean getLoop() { return loop; }
    /**
     * Determine whether the movie includes elevation animation
     * @return true if set to animate
     */
    public boolean getElevationAnimate() { return elevAnim; }
    /**
    * Determine whether the movie includes azimuth animation
    * @return true if set to animate
    */
    public boolean getAzimuthAnimate() { return azimAnim; }
    /**
     * Get the elevation animation start angle
     * @return angle in degrees
     */
    public int getElevationStart() { return elevStart; }
    /**
     * Get the elevation animation increment angle
     * @return increment in degrees
     */
    public int getElevationIncrement() { return elevInc; }
    /**
     * Get the azimuith animation start angle
     * @return angle in degrees
     */
    public int getAzimuthStart() { return azimStart; }
    /**
     * Get the azimuth animation increment angle
     * @return increment in degrees
     */
    public int getAzimuthIncrement() { return azimInc; }
    /**
     * Get the current frame
     * @return frame number (aka slice number)
     */
    public float getCurrentFrame() { return current; }
    /**
     * Get the current azimuth
     * @return azimuth
     */
    public int getCurrentAzimuth() { return currentAzimuth; }
    /**
     * Get the current elevation
     * @return elevation
     */
    public int getCurrentElevation() { return currentElevation; }
    /**
     * Is the movie currently playing
     * @return true if movie is actively playing
     */
    public boolean isPlaying() { return isRunning; }
    /**
     * Get the toolbar instance associated with this movie.
     * @return toolbar instance
     * @see StsMovieActionToolbar(StsModel, boolean)
     */
    public StsToolbar getToolbar() { return toolbar; }

}