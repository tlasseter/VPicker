
/**
 * <p>Title: S2S Development</p>
 * <p>Description: SubVolume Class instantiated by the SubVolume toolbar.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Stuat A. Jackson
 * @version 1.0
 */
package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsSubVolume extends StsMainObject
{
	static final public byte NONE = -1;
	static final public byte SINGLE_SURFACE = 0;
	static final public byte DUAL_SURFACE = 1;
	static final public byte BOX_SET = 2;
	static final public byte WELL_SET = 3;
	static final public byte RESERVOIR_UNIT = 4;

	static final public String[] svTypes = new String[] {"Single Surface", "Dual Surface"};
//		{"Single Surface", "Dual Surface", "Box Set", "Well Set", "Reservoir Block Units"};
	protected boolean isApplied = true;
	protected boolean isInclusive = true;
//	protected StsColor stsColor = new StsColor(Color.RED);
	
    protected byte zDomainOriginal = StsParameters.TD_NONE;
    protected byte zDomainOffset = StsParameters.TD_NONE;
	transient protected StsStringFieldBean zDomainStringBean = new StsStringFieldBean(this, "zDomainString", false, "Original Domain");
    transient protected StsStringFieldBean zDomainOffsetStringBean = new StsStringFieldBean(this, "zDomainOffsetString", true, "Offset Domain");

    static final byte TD_DEPTH = StsProject.TD_DEPTH;
    static final byte TD_TIME = StsProject.TD_TIME;
    static public final float nullValue = StsParameters.nullValue;

	public StsSubVolume()
	{
	}

	public StsSubVolume(boolean persist)
	{
		super(persist);
	}

	public StsSubVolume(boolean persist, String name)
	{
		super(persist);
		setName(name);
	}

	public boolean initialize(StsModel model)
	{
        initializeVisibleFlags();
        zDomainOriginal = model.getProject().getZDomain();
		return true;
	}

	public int getSubVolumeType()
	{
		return type;
	}

/**
	 * SubVolume selected on the Object tree
	 */
	public void treeObjectSelected()
	{
		getSubVolumeClass().selected(this);
		currentModel.win3dDisplayAll();
	}

	static public StsSubVolumeClass getSubVolumeClass()
	{
		return (StsSubVolumeClass) currentModel.getCreateStsClass(StsSubVolume.class);
	}

	static public byte getSubVolumeTypeFromString(String typeString)
	{
		for(int i=0; i<svTypes.length; i++)
			if(typeString.equalsIgnoreCase(svTypes[i]))
				return (byte) i;
		return 0;
	}
	public void setClassVisible()
	{
		StsSubVolumeClass subVolumeClass = getSubVolumeClass();
		subVolumeClass.setIsVisible(true);
	}

	public void setIsInclusive(boolean tf)
	{
        isInclusive = tf;
	}

	public boolean getIsInclusive()
	{
        return isInclusive;
	}

	public byte[] computeBitmapCursor(byte[] plane)
	{return null;
	};

	public void addUnion(byte[] subVolumePlane, int dir, float dirCoordinate, StsRotatedGridBoundingBox boundingBox, byte zDomainData)
	{
		StsException.systemError(this, "addUnion", "Method addUnion() not implemented.");
	}

	public StsGridBoundingBox getGridBoundingBox()
	{
		StsException.systemError(this, "getLoopBoundingBox", "Method getBoundingBox() not implemented.");
		return null;
	}

	public void setIsVisible(boolean isVisible)
	{
		if (this.isVisible == isVisible)return;
		this.isVisible = isVisible;
		resetIsVisible();
		dbFieldChanged("isVisible", isVisible);
		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public void resetIsVisible()
	{
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass) currentModel.getStsClass(getClass());
		subVolumeClass.checkSetClassIsVisible(isVisible);
	}

	/** returns true if changed */
	public boolean setIsVisibleNoDisplay(boolean isVisible)
	{
		if(this.isVisible == isVisible) return false;
		this.isVisible = isVisible;
		return true;
	}

	public void initializeVisibleFlags()
	{
		resetIsVisible();
		resetIsApplied();
	}

	public void setIsApplied(boolean isApplied)
	{
		if (this.isApplied == isApplied) return;
		this.isApplied = isApplied;
		resetIsApplied();
		dbFieldChanged("isApplied", isApplied);	
		currentModel.viewObjectChangedAndRepaint(this, this);
	}

	public void resetIsApplied()
	{
		StsSubVolumeClass subVolumeClass = (StsSubVolumeClass) currentModel.getStsClass(getClass());
		subVolumeClass.checkSetClassIsApplied(isApplied);
	}

	public boolean setIsAppliedNoChange(boolean isApplied)
	{
		if(this.isApplied == isApplied)
            return false;
		this.isApplied = isApplied;
		return true;
	}
	public boolean getIsApplied() { return isApplied; }
    public void display(StsGLPanel panel)  { }
    public byte getZDomain() { return zDomainOriginal; }
    public void setZDomain(byte domain) { zDomainOriginal = domain; }
    public String getZDomainString() { return StsParameters.TD_ALL_STRINGS[zDomainOriginal]; }

    /** the subvolume required is in zDomainData (TIME or DEPTH).  The surface has an original zDomainData: zDomainData.  The surface may have an
     *  offset in a different zDomainData: zDomainOffset.  Compute the correct time or depth value in zDomainData, given this surface and offset.
      * @param surface
     * @param row row index of point
     * @param col col index of point
     * @param offset offset (in offset zDomain)
     * @param velocityModel velocity model used if needed for time-depth conversion
     * @return  time or depth value in zDomainData
     */
    protected float getSurfaceZorT(StsSurface surface, int row, int col, float offset, byte zDomainData, StsSeismicVelocityModel velocityModel)
    {
        float x = surface.getXCoor(row, col);
        float y = surface.getYCoor(row, col);
        float z = surface.getInterpolatedZorT(zDomainData, x, y);

        if(velocityModel == null)
            return z + offset;

        byte zDomainSurfaceOriginal = surface.zDomainOriginal;
        if(zDomainData == zDomainSurfaceOriginal)
        {
            if(zDomainData == zDomainOffset)
                return z + offset;
            else // zDomainData != zDomainOffset && zDomainData == zDomainSurfaceOriginal
            {
                if(zDomainData == TD_TIME) // zDomainOffset == TD_DEPTH
                {
                    z = surface.getZorT(TD_DEPTH, row, col);
                    if(z == nullValue) return z;
                    return (float)velocityModel.getT(x, y, z + offset, 0.0f);
                }
                else // zDomainData == TD_DEPTH && zDomainOffset == TD_TIME
                {
                    float t = surface.getZorT(TD_TIME, row, col);
                    if(t == nullValue) return t;
                    return (float)velocityModel.getZ(x, y, t + offset);
                }
            }
        }
        else // zDomainData != zDomainSurfaceOriginal
        {
            if(zDomainData == TD_TIME) // zDomainSurfaceOriginal == TD_DEPTH
            {
                if(zDomainOffset == TD_DEPTH)
                {
                    z = surface.getZorT(TD_DEPTH, row, col);
                    if(z == nullValue) return z;
                    return (float)velocityModel.getT(x, y, z + offset, 0.0f);
                }
                else // zDomainOffset == TD_TIME
                {
                    return surface.getInterpolatedZorT(TD_TIME, x, y) + offset;
                }
            }
            else // zDomainData == TD_DEPTH && zDomainSurfaceOriginal == TD_TIME
            {
                if(zDomainOffset == TD_DEPTH)
                {
                    return surface.getInterpolatedZorT(TD_DEPTH, x, y) + offset;
                }
                else // zDomainOffset == TD_TIME
                {
                    float t = surface.getInterpolatedZorT(TD_TIME, x, y);
                    if(t == nullValue) return t;
                    return (float)velocityModel.getZ(x, y, t + offset);
                }
            }
        }
    }

    /** the subvolume required is in zDomainData (TIME or DEPTH).  The surface has an original zDomainData: zDomainData.  The surface may have an
     *  offset in a different zDomainData: zDomainOffset.  Compute the correct time or depth value in zDomainData, given this surface and offset.
     * @param surface
     * @param x x-coordinate of point
     * @param y y-coordinate of point
     * @param offset offset (in offset zDomain)
     * @param velocityModel velocity model used if needed for time-depth conversion
     * @return  time or depth value in zDomainData
     */
    protected float getSurfaceZorT(StsSurface surface, float x, float y, float offset, byte zDomainData, StsSeismicVelocityModel velocityModel)
    {
        float z = surface.getInterpolatedZorT(zDomainData, x, y);
        if(z == nullValue) return z;

        if(velocityModel == null)
            return z + offset;;

        byte zDomainSurfaceOriginal = surface.zDomainOriginal;
        if(zDomainData == zDomainSurfaceOriginal)
        {
            if(zDomainData == zDomainOffset)
                return z + offset;
            else // zDomainData != zDomainOffset && zDomainData == zDomainSurfaceOriginal
            {
                if(zDomainData == TD_TIME) // zDomainOffset == TD_DEPTH
                {
                    z = surface.getInterpolatedZorT(TD_DEPTH, x, y);
                    if(z == nullValue) return z;
                    return (float)velocityModel.getT(x, y, z + offset, 0.0f);
                }
                else // zDomainData == TD_DEPTH && zDomainOffset == TD_TIME
                {
                    float t = surface.getInterpolatedZorT(TD_TIME, x, y);
                    if(t == nullValue) return t;
                    return (float)velocityModel.getZ(x, y, t + offset);
                }
            }
        }
        else // zDomainData != zDomainSurfaceOriginal
        {
            if(zDomainData == TD_TIME) // zDomainSurfaceOriginal == TD_DEPTH
            {
                if(zDomainOffset == TD_DEPTH)
                {
                    z = surface.getInterpolatedZorT(TD_DEPTH, x, y);
                    if(z == nullValue) return z;
                    return (float)velocityModel.getT(x, y, z + offset, 0.0f);
                }
                else // zDomainOffset == TD_TIME
                {
                    return surface.getInterpolatedZorT(TD_TIME, x, y) + offset;
                }
            }
            else // zDomainData == TD_DEPTH && zDomainSurfaceOriginal == TD_TIME
            {
                if(zDomainOffset == TD_DEPTH)
                {
                    return surface.getInterpolatedZorT(TD_DEPTH, x, y) + offset;
                }
                else // zDomainOffset == TD_TIME
                {
                    float t = surface.getInterpolatedZorT(TD_TIME, x, y);
                    if(t == nullValue) return t;
                    return (float)velocityModel.getZ(x, y, t + offset);
                }
            }
        }
    }
}
