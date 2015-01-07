package com.Sts.DBTypes;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.StsViewSelectable;
import com.Sts.Interfaces.StsViewable;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import java.io.*;
import java.text.*;
import java.util.*;

public class StsWellPlan extends StsWell implements StsViewable, StsViewSelectable
{
	/** Trajectory type for this well plan */
	protected String trajectoryType = BUILD_HOLD_STRING;
	/** time stamp used in naming this plan */
	protected long time;
	/** platform this well is on */
	protected StsPlatform platform;
	/** From KOP, inclination change in degrees per 100 measured-depth-units.*/
	protected float buildRate = 3.0f;
	/** In drop section, inclination change in degrees per 100 measured-depth-units. */
	protected float dropRate = 3.0f;
	/** In dodge section, inclination change in degrees per 100 measured-depth-units. */
	protected float midPointRate = 3.0f;
	/** for a build, hold, and drop trajectory, this is the distance above the
	 *  target at the end of the dropoff.  From this point the trajectory
	 *  proceeds along the dropoffInclinationAngle into the target. */
	protected float dzAtDropoffEndAboveTarget = 1000.0f;
	/** for a build, hold, and drop trajectory, this is the the angle for
	 *  the final straight-line section into the target in degrees. */
	protected float dropoffInclinationAngle = 45.0f;
	/** accuracy of display and computed wellPlan paths,  used only in arc displays */
	protected float displayInterval = 100.0f;

	/** location of KB in 3D rotated relative coordinates */
	protected StsPoint kbPoint = null;
	/** location of kickoff in 3D rotated relative coordinates */
	protected StsPoint kickoffPoint = null;
	/** location of target in 3D rotated relative coordinates */
	protected StsPoint targetPoint = null;
	/** location of any intermediate points in 3D rotated relative coordinates; includes midPointRate as coordinate #6. */
	protected StsPoint[] midPoints = null;
	/** z coordinate of KellyBushing from sea level; negative if above.  */
	protected float zKB;
	/** z coordinate of kickoff point from sea level; positive down.  */
	protected float zKickoff;
	/** absolute X coordinate of target */
	protected double targetX;
	/** absolute Y coordinate of target */
	protected double targetY;
	/** time coordinate of target */
	protected float targetT;
	/** depth coordinate of target (computed from velocity model) */
	protected float targetZ;
	/** Drilling well this plan may be coming off of (if not null). */
	protected StsWell drillingWell = null;
	/** project coordinates of computed wellPlan */
	protected StsPoint[] wellPlanPoints = null;

	/** Wellplan that this wellplan may have been built from */
	transient StsWellPlan prevWellPlan = null;
	/** absolute X coordinate of midPoint */
	transient protected double midPointX;
	/** absolute Y coordinate of midPoint */
	transient protected double midPointY;
	/** time coordinate of midPoint */
	transient float midPointT;

	/** location of any intermediate points in 3D rotated relative coordinates; includes midPointRate as coordinate #6. */
	transient protected StsPoint[] rotatedMidPoints = null;

//	protected float wellStartZ = -StsParameters.largeFloat;

	/** XYZMT point for last point of actual well (start of new plan) */
	StsPoint wellStartPoint = null;
	/** XYZ vector for last point of actual well (start of new plan) */
	StsPoint wellStartVector = null;
	/** Initial turning rate if from existing well */
	float wellRate = 3.0f;
	/** bottom point of drilling well */
//	transient StsPoint drillingBotPoint;
	/** indicates drilling well is beyond kickoff */
	transient boolean isDrillingBeyondKickoff;

	/** Leg segments for this well plan */
	transient protected Leg[] wellLegsList = null;
	transient StsSurfaceVertex currentVertex = null;
	transient StsProject project;
	transient StsSeismicVelocityModel velocityModel;
    transient SimpleDateFormat fullFormat = new SimpleDateFormat("EEE MMM dd kk:mm:ss z yyyy");

	/** Indicates wellPlan is being built or edited.  While under construction, draw line as solid, othewise dotted. */
	public transient boolean isConstructing = true;

	static public boolean debug = false;

    public static final StsColor defaultWellPlanColor = StsColor.BLUE;

    public static final String STRAIGHT_STRING = "Straight";
	public static final String BUILD_HOLD_STRING = "Build & hold";
	public static final String BUILD_HOLD_DROP_STRING = "Build, hold, & drop";
	public static final String DEEP_KOP_STRING = "Deep kick-off point";
	public static final String HORIZ_WELL_STRING = "Horizontal well";

	public static final String[] trajectoryTypeStrings = new String[]
		{STRAIGHT_STRING, BUILD_HOLD_STRING, BUILD_HOLD_DROP_STRING};

	StsPoint verticalVector = new StsPoint(0.0f, 0.0f, 1.0f, 0.0f, 0.0f);

    public static final float bigFloat = 1.0e10f;

    public StsWellPlan()
	{
	}

	/** Construct wellPlan with persistent == false, and addToModel when successfully constructed. */
	public StsWellPlan(boolean persistent)
	{
		super(getTimeStampName(), persistent, new StsColor(defaultWellPlanColor));
        setZDomainOriginal(StsParameters.TD_TIME_DEPTH);

//		initializePlanPoints();
//		initializeColor();
	}

	/** Construct wellPlan with persistent == false, and addToModel when successfully constructed. */
	public StsWellPlan(String name, boolean persistent)
	{
		super(name, persistent, new StsColor(defaultWellPlanColor));
//		initializePlanPoints();
//		initializeColor();
	}

	public void initializePlanPoints()
	{
		setName("WellPlan-" + getIndex());
//        plannedWells = StsObjectRefList.constructor(2, 2, "plannedWells", this);
		if (currentModel != null)
		{
			project = currentModel.getProject();

			// classInitialize well location to center of cube
			float xCenter = project.getXCenter();
			float yCenter = project.getYCenter();
			double[] xy = project.getAbsoluteXYCoordinates(xCenter, yCenter);
			xOrigin = xy[0];
			yOrigin = xy[1];

			// classInitialize Kelly Bushing to project top depth and
			// kick off depth to 500 below that or at 200 subsea.
			zKB = project.getZorTMin();
			velocityModel = project.getSeismicVelocityModel();
			targetT = project.getTimeMax();
		}
		zKickoff = Math.max(zKB + 500.0f, 200.0f);
		// classInitialize build-rate to 3 deg per hundred
		buildRate = 3.0f;
	}

	public boolean initialize(StsModel model)
	{
		isConstructing = false;
		initializePlan();
		return true;
	}

	public void initializePlan()
	{
		project = currentModel.getProject();
		velocityModel = project.getSeismicVelocityModel();
	}

	public void addTimeStampName()
	{
		time = System.currentTimeMillis();
		Date date = new Date(time);
		String name = date.toString();
		setName(name);
	}

    public long getTimeStamp()
    {
        try
        {
            Date date = fullFormat.parse(getName());
            return date.getTime();
        }
        catch(Exception e)
        {
            return 0l;
        }
    }

	static public String getTimeStampName()
	{
		long time = System.currentTimeMillis();
		Date date = new Date(time);
		String name = date.toString();
		return name;
	}

	public void setPlatform(StsPlatform platform)
	{
		if(platform == null) return;
		this.platform = platform;
		xOrigin = platform.getXOrigin();
		yOrigin = platform.getYOrigin();
	}
	public StsPlatform getPlatform() { return platform; }

	/*
		public boolean delete()
		{
			if (!super.delete())
			{
				return false;
			}
			if(currentPlan != null) currentPlan.delete();
			StsObjectRefList.deleteAll(plannedWells);
			deleteDrawPoints();
			currentPlan = null;
			return true;
		}
	 */
	/** From the specified data, build a set of StsPoints in rotated relative coordinates. */
	public boolean constructPlan()
	{
		initializePoints();
		if (trajectoryType.equals(STRAIGHT_STRING))
			return constructVerticalLeg();
		else if (trajectoryType.equals(BUILD_HOLD_STRING))
			return constructBuildHoldLegs();
		else if (trajectoryType.equals(BUILD_HOLD_DROP_STRING))
			return constructBuildHoldDropLegs();
		else
			return false;
	}

	private void initializePoints()
	{
		clearDrawPoints();

		adjustKBPoint();
		adjustKickoffPoint();
		computeTargetPoint();
		initializeWellLegsList();
    }

	private void adjustKickoffPoint()
	{
		if(kickoffPoint == null)
			kickoffPoint = new StsPoint(5);
		kickoffPoint.v[0] = 0.0f;
		kickoffPoint.v[1] = 0.0f;
		kickoffPoint.v[2] = zKickoff;
		float dXOrigin = (float) (xOrigin - project.getXOrigin()); // well xOrigin offset from project xOrigin
		float dYOrigin = (float) (yOrigin - project.getYOrigin()); // well yOrigin offset from project yOrigin
		rotatePoint(kickoffPoint, dXOrigin, dYOrigin);
		computeTime(kickoffPoint);
	}

	private void computeTargetPoint()
	{
		float dXOrigin = (float) (xOrigin - project.getXOrigin()); // well xOrigin offset from project xOrigin
		float dYOrigin = (float) (yOrigin - project.getYOrigin()); // well yOrigin offset from project yOrigin

		targetPoint = new StsPoint(5);
		targetPoint.setX( (float) targetX); // unrotated relative X
		targetPoint.setY( (float) targetY); // unrotated relative Y
		targetPoint.setT(targetT);
		rotatePoint(targetPoint, dXOrigin, dYOrigin);
		computeDepth(targetPoint);
	}

	public void constructDrillingPlan(boolean isNewPlanSet)
	{
		initializePlan();
		initializePoints();
		checkDrillingPlan();
		if(!isNewPlanSet) updateDrillingPlan();
	}

	public void updateDrillingPlan()
	{
		if (trajectoryType.equals(BUILD_HOLD_STRING))
		{
			if(!constructDrillingBuildHoldLegs()) return;
		}
		else if (trajectoryType.equals(BUILD_HOLD_DROP_STRING))
		{
			if(!constructDrillingBuildHoldDropLegs()) return;
		}
	}

	public boolean computePoints()
	{
		int nIncrements = 0;
		int nWellLegs = wellLegsList.length;
		StsPoint[][] wellLegPoints = new StsPoint[nWellLegs][];
		for (int n = 0; n < nWellLegs; n++)
		{
			wellLegPoints[n] = wellLegsList[n].getPoints();
			nIncrements += wellLegPoints[n].length - 1;
		}
		wellPlanPoints = new StsPoint[nIncrements + 1];
		int p = 0;
		for (int n = 0; n < nWellLegs; n++)
		{
			int nPoints = 0;
			if (n == nWellLegs - 1)
			{
				nPoints = wellLegPoints[n].length;
			}
			else
			{
				nPoints = wellLegPoints[n].length - 1;
			}

			for (int i = 0; i < nPoints; i++)
			{
				wellPlanPoints[p++] = wellLegPoints[n][i];
			}
		}
		StsLine.addMDepthToPoints(wellPlanPoints);
		return true;
	}
	public StsPoint[] getExportPoints()
	{
        return this.wellPlanPoints;
    }
/*
	public StsPoint[] getExportPoints()
	{
        velocityModel = project.getSeismicVelocityModel();

        int nIncrements = 0;
		if(wellLegsList == null)
		{
			if(!constructPlan())return null;
			insertMidPoints();
		}
		int nWellLegs = wellLegsList.length;
		StsPoint[][] wellLegExportPoints = new StsPoint[nWellLegs][];
		for (int n = 0; n < nWellLegs; n++)
		{
			wellLegExportPoints[n] = wellLegsList[n].getExportPoints();
			nIncrements += wellLegExportPoints[n].length - 1;
		}
		StsPoint[] exportPoints = new StsPoint[nIncrements + 1];
		int p = 0;
		for (int n = 0; n < nWellLegs; n++)
		{
			int nPoints = 0;
			if (n == nWellLegs - 1)
			{
				nPoints = wellLegExportPoints[n].length;
			}
			else
			{
				nPoints = wellLegExportPoints[n].length - 1;
			}

			for (int i = 0; i < nPoints; i++)
			{
				exportPoints[p++] = wellLegExportPoints[n][i];
			}
		}
		StsLine.addMDepthToPoints(exportPoints);
		return exportPoints;
	}
*/
	public StsPoint[] getRotatedPoints() { return wellPlanPoints; }

	public void clearPoints()
	{
		wellPlanPoints = null;
	}

	public boolean constructVerticalLeg()
	{
		try
		{
			StartStraightLeg verticalLeg = new StartStraightLeg("Vertical leg", kbPoint, targetPoint); // defines a vertical kb to target time leg
			if(!doLegs()) return false;
			computePoints();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.initializeVerticalLeg() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private boolean constructBuildHoldLegs()
	{
		try
		{
			float radius = (float) (18000 / (Math.PI * buildRate));
			StartStraightLeg kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
			StartArcLeg kickoffArcLeg = new StartArcLeg("Kickoff arc", radius);
			EndStraightLeg holdLeg = new EndStraightLeg("Target leg", targetPoint);
			if(!doLegs()) return false;
			computePoints();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.initializeBuildHoldLegs() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private boolean constructDrillingBuildHoldLegs()
	{
		try
		{
			if(!isDrillingBeyondKickoff)
			{
				StartStraightLeg kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
				float radius = (float) (18000 / (Math.PI * buildRate));
				StartArcLeg kickoffArcLeg = new StartArcLeg("Kickoff arc", radius);
			}
			else
			{
				float wellStartRadius = (float) (18000 / (Math.PI * wellRate));
				StartWellArcLeg startWellLeg = new StartWellArcLeg("Start well Leg", wellStartRadius, wellStartPoint, wellStartVector);
            }
			EndStraightLeg holdLeg = new EndStraightLeg("Target leg", targetPoint);
            // hacque: figure out what to do!
//            holdLeg.initializeInVector(wellStartPoint);
			if (!doLegs())return false;
			computePoints();
			return insertMidPoints();
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.initializeBuildHoldLegs() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private boolean insertMidPoints()
	{
		if(midPoints == null) return true;
		int nPoints = midPoints.length;
		StsPoint[] newMidPoints = midPoints;
		midPoints = null;
		for (int n = 0; n < nPoints; n++)
			computeMidPoint(newMidPoints[n]);
		return true;
	}

	private boolean constructDrillingBuildHoldDropLegs()
	{
		try
		{
			if(!isDrillingBeyondKickoff)
			{
				StartStraightLeg kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
				float radius = (float) (18000 / (Math.PI * buildRate));
				StartArcLeg kickoffArcLeg = new StartArcLeg("Kickoff arc", radius);
			}
			else
			{
				float wellStartRadius = (float) (18000 / (Math.PI * wellRate));
				StartWellArcLeg startWellLeg = new StartWellArcLeg("Start well Leg", wellStartRadius, wellStartPoint, wellStartVector);
			}
			MidStraightLeg holdLeg = new MidStraightLeg("Hold leg");
			float dropRadius = (float) (18000 / (Math.PI * dropRate));
			EndArcLeg dropArcLeg = new EndArcLeg("Drop arc", dropRadius);
			InclinedEndLeg targetLeg = new InclinedEndLeg("Inclined target leg", dropoffInclinationAngle, dzAtDropoffEndAboveTarget, targetPoint);
			if(!doLegs()) return false;
			computePoints();
			return insertMidPoints();
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.initializeBuildHoldLegs() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private void initializeWellLegsList()
	{
		wellLegsList = new Leg[0];
	}

	/** Create a newArcLeg which splits the existing holdLeg into prevLeg and newHoldLeg.
	 *  Relink holdLeg, newArcLeg, newHoldLeg, and nextLeg (if it exists).
	 */
	private Leg[] insertArcLegInStraightLeg(String name, float radius, StsPoint midTargetPoint, Leg holdLeg)
	{
		MidArcLeg newArcLeg = new MidArcLeg(name, radius, midTargetPoint);
		holdLeg.legName = "Hold leg after " + name;
		wellLegsListInsertBefore(newArcLeg, holdLeg);
		Leg newHoldLeg = constructStraightLegInsertBefore("Hold leg before " + name, newArcLeg);
		return new Leg[]
			{newHoldLeg, newArcLeg};
	}

	/** we are assuming that two legs, one arc one straight are to be deleted. */
	private void deleteLegs(Leg[] deleteLegs)
	{
		if (deleteLegs.length != 2)
		{
			StsException.systemError("StsWellPlanSet.deleteLegs() failed. Number of legs is not 2.");
			return;
		}
		if (debug)
		{
			System.out.println("    Deleted legs " + deleteLegs[0].legName + " and " + deleteLegs[1].legName);
		}
		deleteLeg(deleteLegs[0]);
		deleteLeg(deleteLegs[1]);
	}

	private void deleteLeg(Leg leg)
	{
		wellLegsList = (Leg[]) StsMath.arrayDeleteElement(wellLegsList, leg);
	}

	private boolean constructBuildHoldDropLegs()
	{
		try
		{
			float buildRadius = (float) (18000 / (Math.PI * buildRate));
			float dropRadius = (float) (18000 / (Math.PI * dropRate));
			StartStraightLeg kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
			StartArcLeg kickoffArcLeg = new StartArcLeg("Kickoff arc", buildRadius);
			MidStraightLeg holdLeg = new MidStraightLeg("Hold leg");
			EndArcLeg dropArcLeg = new EndArcLeg("Drop arc", dropRadius);
			InclinedEndLeg targetLeg = new InclinedEndLeg("Inclined target leg", dropoffInclinationAngle, dzAtDropoffEndAboveTarget, targetPoint);
			if(!doLegs()) return false;
			computePoints();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.initializeBuildHoldDropLegs() failed.", e, StsException.WARNING);
			return false;
		}
	}

	private boolean doLegs()
	{
		initializeLegs();
		setupLegs();
		if (debug) debugLegs("\nAFTER SETUP\n");
		computeLegs();
        if (debug) debugLegs("\nAFTER COMPUTE\n");
        boolean adjustedOK = adjustLegs();
        if (debug) debugLegs("\nAFTER ADJUST\n");
        if(adjustedOK)
            completeLegs();
		return adjustedOK;
	}

	private void initializeLegs()
	{

		int nLegs = wellLegsList.length;
		Leg prevLeg = null;
		for (int n = 0; n < nLegs; n++)
		{
			Leg leg = wellLegsList[n];
			leg.initialize();
			leg.prevLeg = prevLeg;
			if (prevLeg != null)
			{
				prevLeg.nextLeg = leg;
			}
			prevLeg = leg;
		}
	}

	private void setupLegs()
	{
		int nLegs = wellLegsList.length;
		for (int n = 0; n < nLegs; n++)
		{
			wellLegsList[n].linkLegs();
		}
	}

	private void computeLegs()
	{
		int nLegs = wellLegsList.length;
		wellLegsList[0].compute();
		wellLegsList[nLegs - 1].compute();
		for (int n = 1; n < nLegs - 1; n++)
		{
//        for(int n = 0; n < nLegs; n++)
			wellLegsList[n].compute();
		}
	}

	private boolean adjustLegs()
	{
		int nIters = 0;

		int nLegs = wellLegsList.length;
		boolean adjusting = true;
		while (adjusting)
		{
			adjusting = false;
			for (int n = 0; n < nLegs; n++)
			{
				if (wellLegsList[n].adjust())
				{
					adjusting = true;
				}
			}
			nIters++;
			if (nIters > 100)
			{
                StsException.systemError(this, "adjustLegs", "Iterations exceeded 100.");
                for (int n = 0; n < nLegs; n++)
                {
                    if (!wellLegsList[n].converged)
                    {
                        System.out.println(wellLegsList[n].legName + " failed to converged");
                    }
                }
                return false;
			}
		}
		return true;
	}

	private void completeLegs()
	{
		int nLegs = wellLegsList.length;
		for (int n = 0; n < nLegs; n++)
		{
			wellLegsList[n].complete();
		}
	}

	private void debugLegs(String string)
	{
        System.out.println(string);
        int nLegs = wellLegsList.length;
		for (int n = 0; n < nLegs; n++)
		{
			wellLegsList[n].debug();
		}
	}

	private void rotatePoint(StsPoint point, float dXOrigin, float dYOrigin)
	{
		float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY(dXOrigin + point.v[0], dYOrigin + point.v[1]);
		point.v[0] = xy[0];
		point.v[1] = xy[1];
	}

	public StsPoint computeTimeDepthPoint(StsPoint point)
	{
		try
		{
			StsPoint timeDepthPoint = new StsPoint(5);
			float x = point.getX();
			float y = point.getY();
			timeDepthPoint.setX(x);
			timeDepthPoint.setY(y);
			StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
			if (isDepth)
			{
				float z = point.getZ();
				timeDepthPoint.setZ(z);
				if (velocityModel == null)return timeDepthPoint;
				float t = (float) velocityModel.getT(point.v);
				timeDepthPoint.setT(t);
			}
			else
			{
				float t = point.getZ();
				timeDepthPoint.setT(t);
				if (velocityModel == null)return timeDepthPoint;
				float z = (float) velocityModel.getZ(x, y, t);
				timeDepthPoint.setZ(z);
			}
			return timeDepthPoint;
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.computeTimeDepthPoint() failed.", e, StsException.WARNING);
			return null;
		}
	}

	public void computeDepth(StsPoint point)
	{
		if (point.getLength() < 5)
		{
			StsException.systemError("StsWellPlanSet.computeDepth() failed.  Received point argument with length < 5.");
			return;
		}
		float x = point.getX();
		float y = point.getY();
		float t = point.getT();
		float z = (float) velocityModel.getZ(x, y, t);
		point.setZ(z);
	}

	public void computeTime(StsPoint point)
	{
		try
		{
			if (velocityModel == null)
			{
				return;
			}

			if (point.getLength() < 5)
			{
				StsException.systemError("StsWellPlanSet.computeDepth() failed.  Received point argument with length < 5.");
				return;
			}
			float t = (float) velocityModel.getT(point.v);
			point.setT(t);
		}
		catch (Exception e)
		{
			StsException.outputException("StsWellPlanSet.computeTime() failed for point " + point.toString(), e, StsException.WARNING);
		}
	}

	/*
		public StsWellPlan copyLastPlannedWell()
		{
			try
			{
				if (plannedWells == null)
				{
					return null;
				}
				StsWellPlan oldPlan = (StsWellPlan) plannedWells.getLast();
				StsWellPlan newPlan = new StsWellPlan();
				//           StsToolkit.copySubToSuperclass(oldPlan, newPlan, StsWellPlan.class, StsMainObject.class, true);

				StsPoint[] oldPoints = oldPlan.getPoints();
				int nPoints = oldPoints.length;
				StsPoint[] newPoints = new StsPoint[nPoints];
				for (int n = 0; n < nPoints; n++)
				{
					newPoints[n] = new StsPoint(oldPoints[n]);
				}
				newPlan.setPoints(newPoints);
				addPlannedWell(newPlan);
				return newPlan;
			}
			catch (Exception e)
			{
				StsException.outputException("StsWellPlanSet.copyLastPlannedWell() failed.",
											 e, StsException.WARNING);
return null;
			}
		}
	 */
	public boolean resetDrillingWell()
	{
		drillingWell = (StsWell) currentModel.getObjectWithName(StsWell.class, drillingWell.getName());
		if (drillingWell != null)
		{
			return true;
		}
		return false;
	}

	/*
		public void addAbsolutePoint(double x, double y, double z, double m, double t, float rate)
		{
			StsSurfaceVertex vertex = currentPlan.addAbsolutePoint(x, y, z, m, t);
//        currentPlan.adjustRates(currentPlan.getVertexIndex(vertex), rate, currentPlan.INSERT);
		}

		public void addAbsolutePoint(double x, double y, double z, float rate)
		{
			StsSurfaceVertex vertex = currentPlan.addAbsolutePoint(x, y, z);
			if (vertex != null)
			{
//            currentPlan.adjustRates(currentPlan.getVertexIndex(vertex), rate, currentPlan.INSERT);
			}
		}
	 */
	public void setTargetPoint(StsPoint point)
	{
		setTargetX(point.getX());
		setTargetY(point.getY());
		setTargetT(point.getT());
		computeTargetPoint();
	}

	public void setTargetX(double x)
	{
		targetX = x;
	}

	public void setTargetY(double y)
	{
		targetY = y;
	}

	public void setTargetT(float t)
	{
		targetT = t;
	}

	public double getTargetX()
	{
		return targetX;
	}

	public double getTargetY()
	{
		return targetY;
	}

	public float getTargetT()
	{
		return (float) targetT;
	}

	public void setMidPoint(StsPoint point)
	{
		setMidPointX(point.getX());
		setMidPointY(point.getY());
		setMidPointT(point.getT());
	}

	public void setMidPointX(double x)
	{
		midPointX = x;
	}

	public void setMidPointY(double y)
	{
		midPointY = y;
	}

	public void setMidPointT(float t)
	{
		midPointT = t;
	}

	public void setMidPointRate(float rate)
	{
		midPointRate = rate;
	}

	public double getMidPointX()
	{
		return midPointX;
	}

	public double getMidPointY()
	{
		return midPointY;
	}

	public float getMidPointT()
	{
		return midPointT;
	}

	public float getMidPointRate()
	{
		return midPointRate;
	}

	public float getWellX()
	{
		if (wellStartPoint == null)return bigFloat;
		else return wellStartPoint.getX();
	}

	public void setWellX(float x)
	{
		wellStartPoint.setX(x);
	}

	public float getWellY()
	{
		if (wellStartPoint == null)return bigFloat;
		else return wellStartPoint.getY();
	}

	public void setWellY(float y)
	{
		wellStartPoint.setX(y);
	}

	public float getWellZ()
	{
		if (wellStartPoint == null)return bigFloat;
		else return wellStartPoint.getZ();
	}

	public void setWellZ(float z)
	{
		wellStartPoint.setX(z);
	}

	public float getWellT()
	{
		if (wellStartPoint == null)return bigFloat;
		else return wellStartPoint.getT();
	}

	public void setWellT(float t)
	{
		wellStartPoint.setX(t);
	}

	public float getWellRate()
	{
		return wellRate;
	}

	public void setWellRate(float rate)
	{
		wellRate = rate;
	}

	public void deleteLineVertex(StsSurfaceVertex deleteVertex)
	{
		int vIdx = getVertexIndex(deleteVertex);
		deleteLineVertex(deleteVertex);
		checkChangeCurrentVertex();
//        currentPlan.adjustRates(vIdx, 0.0f, currentPlan.DELETE);
		computePoints();
	}

	public StsWell getDrillingWell()
	{
		return drillingWell;
	}

	public void setDrillingWell(StsWell well)
	{
		drillingWell = well;
		wellStartPoint = well.getBotRotatedPoint();
		wellStartVector = well.getBotVector();
//		xOrigin = well.getXOrigin();
//		yOrigin = well.getYOrigin();
		return;
	}

	public void setPrevWellPlan(StsWellPlan prevPlan)
	{
		prevWellPlan = prevPlan;
	}

	public boolean checkDrillingPlan()
	{
		String units = project.getXyUnitString();
//		drillingBotPoint = drillingWell.getBotPoint();
		if (!checkKBMatch(units))return false;
		if (!checkKickoffMatch(units))return false;
		return true;
	}

	private boolean checkKBMatch(String units)
	{
		StsPoint wellKbPoint = this.drillingWell.getTopPoint();
    /*
		float error = Math.abs(wellKbPoint.getZ() - kbPoint.getZ());
		if (error > 10.0f)
		{
			boolean ok = StsYesNoDialog.questionValue(currentModel.win3d, "Well kb and plan kb differ by " + error + " " + units +
												  ".\n Is this OK? If not, process will be terminated.");
			return ok;
		}
     */
		return true;
	}

	private boolean checkKickoffMatch(String units)
	{
		boolean ok = true;
		float wellMDepth = wellStartPoint.getM();
		float wellDepth = wellStartPoint.getZ();
		float estKickoffMDepth = zKickoff - zKB;
		if(wellMDepth <= estKickoffMDepth)
		{
			if(wellDepth > zKickoff)
			{
				ok = StsYesNoDialog.questionValue(currentModel.win3d, "Depth and measured depth are in error. Well mDepth " + wellMDepth + " is above plan kickoff " + estKickoffMDepth +
											  " but well depth  " + wellDepth + " is greater than plan kickoff depth " + zKickoff +
											  ".\n Is this OK? If not, process will be terminated.");
			}
			isDrillingBeyondKickoff = false;
			return ok;
		}
		float zWell = drillingWell.getDepthFromMDepth(estKickoffMDepth);
		float error = Math.abs(zWell - zKickoff);
		if (error > 10.0f)
		{
            /*
			ok = StsYesNoDialog.questionValue(currentModel.win3d, "Well kickoff and plan kickoff are different by  " + error + " " + units +
												  ".\n Is this OK? If not, process will be terminated.");
            */
		}
		isDrillingBeyondKickoff = true;
		return true;
	}

/*
    public void setStartZ(float z)
    {
        wellStartZ = z;
    }
*/
/*
    public void addPlannedWell(StsWellPlan plannedWell)
    {
        plannedWells.add(plannedWell);
        plannedWell.addTimeStampName();
        plannedWell.setZDomainSupported(StsLine.TD_TIME_DEPTH);
        this.currentPlan = plannedWell;
        plansListBean.addItem(plannedWell);
        StsObjectPanel objectPanel = getObjectPanel();
        if (objectPanel.getPanelObject() == this)
        {
            objectPanel.refreshProperties();
        }
    }

    public Object[] getWellPlans()
    {
        return plannedWells.getTrimmedList();
    }
*/
/*
    public StsFieldBean[] getDisplayFields()
    {
        setPlansListBeanList();
        return displayFields;
    }
 */
/*
    private void setPlansListBeanList()
    {
        plansListBean.removeAll();
        int nWellPlans = plannedWells.getSize();
        for (int n = 0; n < nWellPlans; n++)
        {
            plansListBean.addItem(plannedWells.getElement(n));
        }
        plansListBean.setValueObject(currentPlan);
    }
*/
    public StsObjectRefList getVertices()
    {
        return getLineVertices();
    }
/*
    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }
*/
    /** returns true if changed */
    public boolean setIsVisibleNoDisplay(boolean isVisible)
    {
        if (this.isVisible == isVisible)
        {
            return false;
        }
        this.isVisible = isVisible;
        return true;
    }
/*
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
        currentModel.getStsClass(StsWellPlan.class).selected(this);
    }

    public boolean export()
    {
        if (this.currentPlan == null)
        {
            return false;
        }
        return currentPlan.export(StsParameters.TD_TIME_DEPTH_STRING);
    }
*/
    public void display(StsGLPanel3d glPanel3d, boolean displayName, String drawLabelString, float labelInterval)
	{
        labelFormat = new DecimalFormat(getWellClass().getLabelFormatAsString());
		if (!currentModel.getProject().canDisplayZDomain(zDomainSupported))
		{
			return;
		}
		if (glPanel3d == null)
		{
			return;
		}

	    if (debug && wellLegsList != null)
		{
			GL gl = glPanel3d.getGL();
			int zIndex = StsPoint.getVerticalIndex();
			int nLegs = wellLegsList.length;
			for (int n = 0; n < nLegs; n++)
			{
				wellLegsList[n].display(gl, stsColor, false, zIndex, this);
			}
		}
		if (isConstructing)
		{
			if (kickoffPoint != null)
			{
				StsGLDraw.drawPoint(kickoffPoint.getXYZorT(), stsColor, glPanel3d, 8);
			}
			if (targetPoint != null)
			{
				StsGLDraw.drawPoint(targetPoint.getXYZorT(), stsColor, glPanel3d, 8);
			}
			if (rotatedMidPoints != null)
			{
				for (int n = 0; n < rotatedMidPoints.length; n++)
				{
					StsGLDraw.drawPoint(rotatedMidPoints[n].getXYZorT(), stsColor, glPanel3d, 8);
				}
			}
		}

		GL gl = glPanel3d.getGL();

		if (isDrawingCurtain)
        {
            glPanel3d.setViewShift(gl, 4.0f);
        }
		if (displayName)
			display(glPanel3d, isConstructing, getName(), wellPlanPoints, !isConstructing);
		else
			display(glPanel3d, isConstructing, null, wellPlanPoints, !isConstructing);

		if (isDrawingCurtain)
        {
            glPanel3d.resetViewShift(gl);
        }
		if (!drawLabelString.equals(NO_LABEL) && labelInterval >= 0.0f)
		{
			StsPoint point = null;
			float md = 0.0f;
			String label = null;
			int nLabels = (int) (getMaxMDepth() / labelInterval);

			if (isDrawingCurtain)
			{
				StsColor.BLACK.setGLColor(gl);
				glPanel3d.setViewShift(gl, 10.0f);
			}
			else
			{
				stsColor.setGLColor(gl);
				glPanel3d.setViewShift(gl, 1.0f);
			}

			GLBitmapFont font = GLHelvetica12BitmapFont.getInstance(gl);
			int numChars = font.getNumChars();
			for (int i = 0; i < nLabels; i++, md += labelInterval)
			{
				point = getPointAtMDepth( (float) (i * labelInterval), true);
				float[] xyz = point.getXYZorT();
				if ( (md % (5.0f * labelInterval)) != 0.0f)
				{
					StsGLDraw.drawPoint(xyz, null, glPanel3d, 5, 1, 0.0f);
				}
				else
				{
					StsGLDraw.drawPoint(xyz, null, glPanel3d, 10, 2, 0.0f);
					float value = 0.0f;
					if (drawLabelString.equals(MDEPTH_LABEL))
					{
						value = md;
//                       label = Float.toString(md);
					}
					else if (drawLabelString.equals(DEPTH_LABEL))
					{
						value = point.getZ();
//                        label = Float.toString(point.getZ());
					}
					else if (drawLabelString.equals(TIME_LABEL))
					{
						value = point.getT();
//                       label = Float.toString(point.getT());
					}
					label = labelFormat.format(value);
					StsGLDraw.fontOutput(gl, xyz, label, font);
				}
			}
			glPanel3d.resetViewShift(gl);
		}
        glPanel3d.setViewShift(gl, -4.0f);
        displaySeismicCurtain(glPanel3d);
        glPanel3d.resetViewShift(gl);
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayNames)
    {
        if (debug && wellLegsList != null)
        {
            GL gl = glPanel3d.getGL();
            int zIndex = StsPoint.getVerticalIndex();
            int nLegs = wellLegsList.length;
            for (int n = 0; n < nLegs; n++)
            {
                wellLegsList[n].display(gl, stsColor, true, zIndex, this);
            }
        }
		if(isConstructing)
		{
			if (kickoffPoint != null)
			{
				StsGLDraw.drawPoint(kickoffPoint.getXYZorT(), stsColor, glPanel3d, 8);
			}
			if (targetPoint != null)
			{
				StsGLDraw.drawPoint(targetPoint.getXYZorT(), stsColor, glPanel3d, 8);
			}
			if (rotatedMidPoints != null)
			{
				for (int n = 0; n < rotatedMidPoints.length; n++)
				{
					StsGLDraw.drawPoint(rotatedMidPoints[n].getXYZorT(), stsColor, glPanel3d, 8);
				}
			}
		}

        if(wellPlanPoints != null)
        {
            if(displayNames)
                display(glPanel3d, true, getName(), wellPlanPoints, !isConstructing);
            else
                display(glPanel3d, true, null, wellPlanPoints, !isConstructing);
        }
    }


    public void pick(GL gl, StsGLPanel glPanel)
    {
        int zIndex = StsPoint.getVerticalIndex();
        StsGLDraw.pickLine(gl, stsColor, highlighted, wellPlanPoints, zIndex);
    }

    private void displayLegs(GL gl)
    {
        if (wellLegsList == null) return;
        int nLegs = wellLegsList.length;
        int zIndex = StsPoint.getVerticalIndex();
        for (int n = 0; n < nLegs; n++)
        {
            wellLegsList[n].display(gl, stsColor, false, zIndex, this);
        }
    }

    public void setTrajectoryType(String type)
    {
        if (trajectoryType == type)
        {
            return;
        }
        if (wellLegsList != null)
        {
            if (trajectoryType == BUILD_HOLD_STRING && type == BUILD_HOLD_DROP_STRING)
            {
                removeBuildHoldLeg();
 //               addBuildHoldDropLeg();
            }
            else if (trajectoryType == BUILD_HOLD_DROP_STRING && type == BUILD_HOLD_STRING)
            {
                removeBuildHoldDropLeg();
//                addBuildHoldLeg();
            }
        }
        trajectoryType = type;
    }

    private void removeBuildHoldLeg()
    {
        for (int n = 0; n < wellLegsList.length; n++)
        {
            Leg leg = wellLegsList[n];
            if (leg instanceof EndStraightLeg)
            {
                deleteLeg(leg);
            }
        }
    }

    private void addBuildHoldDropLeg()
    {
        InclinedEndLeg targetLeg = new InclinedEndLeg("Inclined target leg", dropoffInclinationAngle, dzAtDropoffEndAboveTarget, targetPoint);
        addLeg(targetLeg);
    }

    private void removeBuildHoldDropLeg()
    {
        for (int n = 0; n < wellLegsList.length; n++)
        {
            Leg leg = wellLegsList[n];
            if (leg instanceof InclinedEndLeg)
            {
                deleteLeg(leg);
            }
        }
    }

    private void addBuildHoldLeg()
    {
        EndStraightLeg targetLeg = new EndStraightLeg("Target leg", targetPoint);
        addLeg(targetLeg);
    }

    private void addLeg(Leg leg)
    {
        wellLegsList = (Leg[]) StsMath.arrayAddElement(wellLegsList, leg);
    }

    public String getTrajectoryType()
    {
        return trajectoryType;
    }

    public void setXOrigin(double x)
    {
        xOrigin = x;
        adjustKBPoint();
    }

    public void setYOrigin(double y)
    {
        yOrigin = y;
        adjustKBPoint();
    }

    public void setZKB(float zKB)
    {
        this.zKB = -zKB;
        adjustKBPoint();
    }

    public void adjustKBPoint()
    {
        if (kbPoint == null)
        {
            kbPoint = new StsPoint(5);
        }
		if(project == null) project = currentModel.getProject();
        float dXOrigin = (float) (xOrigin - project.getXOrigin()); // well xOrigin offset from project xOrigin
        float dYOrigin = (float) (yOrigin - project.getYOrigin()); // well yOrigin offset from project yOrigin

        kbPoint.v[0] = 0.0f;
        kbPoint.v[1] = 0.0f;
        kbPoint.v[2] = zKB;
        rotatePoint(kbPoint, dXOrigin, dYOrigin);
        computeTime(kbPoint);
    }

    public void setZKickoff(float zKickoff)
    {
        this.zKickoff = zKickoff;
    }

    public void setBuildRate(float rate)
    {
        buildRate = rate;
    }

    public void setDropRate(float rate)
    {
        dropRate = rate;
    }

    public void setDzAboveTarget(float dz)
    {
        this.dzAtDropoffEndAboveTarget = dz;
    }

    public void setDropAngle(float angle)
    {
        dropoffInclinationAngle = angle;
    }

    public double getXOrigin()
    {
        return xOrigin;
    }

    public double getYOrigin()
    {
        return yOrigin;
    }

    public float getZKB()
    {
        return zKB;
    }

    public float getBuildRate()
    {
        return buildRate;
    }

    public float getDropRate()
    {
        return dropRate;
    }

    public float getDzAboveTarget()
    {
        return dzAtDropoffEndAboveTarget;
    }

    public float getDropAngle()
    {
        return dropoffInclinationAngle;
    }

    public void addTargetPoint(StsPoint timeDepthPoint)
    {
        resetToUnrotatedRelativeXY(timeDepthPoint); // vertex is in unrotated relative coordinates
        setTargetPoint(timeDepthPoint);
		if(drillingWell != null)
			updateDrillingPlan();
		else
			constructPlan();
    }

	public void computeMidPoint(StsPoint timeDepthPoint)
	{
		StsPoint projectPoint = new StsPoint(timeDepthPoint); // save a copy
        Leg nearestLeg = null;
        float nearestDistance = bigFloat;
        float f;
		if(wellLegsList == null)
		{
			StsException.systemError("StsWellPlan.addMidPoint() failed: wellLegsList is null.");
			return;
		}
        int nLegs = wellLegsList.length;
        for (int n = 0; n < nLegs; n++)
        {
            Leg wellLeg = wellLegsList[n];
            if (wellLeg instanceof MidStraightLeg || wellLeg instanceof EndStraightLeg)
            {
                DistanceToLine distanceToLine = new DistanceToLine(wellLeg.inPoint, wellLeg.outPoint, timeDepthPoint);
                if (!distanceToLine.isPointNear())
                {
                    continue;
                }
                if (distanceToLine.distance < nearestDistance)
                {
                    nearestDistance = distanceToLine.distance;
                    nearestLeg = wellLeg;
                }
            }
        }
        if (nearestLeg == null)
        {
            /*
            new StsMessage(currentModel.win3d, StsMessage.ERROR,
                           "Added point not on any plane normal to existing mid straight leg or end straight leg.\n" +
                           "Note that the point will not be inserted in an existing inclined target leg.\n" +
                           " Point will be ignored.");
            */
            return;
        }

        float dodgeRadius = (float) (18000 / (Math.PI * midPointRate));
        resetToUnrotatedRelativeXY(timeDepthPoint); // vertex is in unrotated relative coordinates
//        setMidPoint(timeDepthPoint);
        float dXOrigin = (float) (xOrigin - project.getXOrigin()); // well xOrigin offset from project xOrigin
        float dYOrigin = (float) (yOrigin - project.getYOrigin()); // well yOrigin offset from project yOrigin
        rotatePoint(timeDepthPoint, dXOrigin, dYOrigin);
        computeDepth(timeDepthPoint);
        Leg[] addedLegs = insertArcLegInStraightLeg("Dodge arc", dodgeRadius, timeDepthPoint, nearestLeg);
        if (doLegs())
        {
            addMidDrawPoint(timeDepthPoint);
			addMidPoint(projectPoint);
			computePoints();
        }
        else
        {
            deleteLegs(addedLegs);
//                clearDrawPoints();
//                clearPoints();
            doLegs(); // redo the legs construction after bad leg has been deleted
			computePoints();
            new StsMessage(currentModel.win3d, StsMessage.ERROR, "Could not find a path solution with this point added.  Will ignore.");

        }
//            completeWellPlan();
    }

    /** add a midPoint to the list.  The point size is 6: 5 coordinates (x,y,z,md,t) and turn rate. */
    private void addMidDrawPoint(StsPoint midPoint)
    {
        StsPoint midDrawPoint = new StsPoint(6);
        midDrawPoint.copyFrom(midPoint);
        midDrawPoint.v[5] = midPointRate;
        if (rotatedMidPoints == null)
        {
            rotatedMidPoints = new StsPoint[] { midDrawPoint };
        }
        else
        {
            rotatedMidPoints = (StsPoint[]) StsMath.arrayAddElement(rotatedMidPoints, midDrawPoint);
        }
    }

	private void addMidPoint(StsPoint inputMidPoint)
	{
		StsPoint midPoint = new StsPoint(6);
		midPoint.copyFrom(inputMidPoint);
		midPoint.v[5] = midPointRate;
		if (rotatedMidPoints == null)
		{
			midPoints = new StsPoint[] { midPoint };
		}
		else
		{
			midPoints = (StsPoint[]) StsMath.arrayAddElement(rotatedMidPoints, midPoint);
		}
    }

    class DistanceToLine
    {
        float f;
        float distance;
        DistanceToLine(StsPoint linePoint1, StsPoint linePoint2, StsPoint point)
        {
            StsPoint line = new StsPoint(3);
            line.subPoints(linePoint2, linePoint1);
            StsPoint vector = new StsPoint(3);
            vector.subPoints(point, linePoint1);
            float dot = line.dot(vector);
            float lineLengthSq = line.lengthSq();
            f = dot / lineLengthSq;
            float vectorLengthSq = vector.lengthSq();
            distance = (float) Math.sqrt(vectorLengthSq - f * dot);
        }

        boolean isPointNear()
        {
            return f >= 0.0f && f <= 1.0f;
        }
    }

    public void wizardFinish()
    {
//        savePointsToDB();
        clearDrawPoints();
//		addToModel();
		isConstructing = false;
    }

    /** For now we don't need to call this as plannedWellSet is constructed and saved in this same session.
     *  Later, if we provide editing capability after plannedWellSet is constructed, then call this method
     *  to save new plan points.
     */

    private void savePointsToDB()
    {
        dbFieldChanged("kbPoint", kbPoint);
        dbFieldChanged("kickoffPoint", kickoffPoint);
        dbFieldChanged("targetPoint", targetPoint);
        dbFieldChanged("midDrawPoints", midPoints);
    }

    public void clearDrawPoints()
    {
        wellLegsList = null;
    }

    public void deleteDrawPoints()
    {
        kbPoint = null;
        kickoffPoint = null;
        targetPoint = null;
        midPoints = null;
        wellLegsList = null;
    }

    /** Given a point in the relative rotated coordinate system (like a picked point),
     *  set the xy coordinates to the unrotated relative coordinates of this StsLine.
     *  These coordinates are unrotated xy coordinate offsets from the origin of the well
     *  which is typically the top point (Kelly bushing) xy.
     */
    public void resetToUnrotatedRelativeXY(StsPoint point)
    {
        double[] xy = currentModel.getProject().getAbsoluteXYCoordinates(point);
        point.v[0] = (float) (xy[0] - xOrigin);
        point.v[1] = (float) (xy[1] - yOrigin);
    }

    public void setDisplayInterval(float displayInterval)
    {
        this.displayInterval = displayInterval;
    }

    public float getDisplayInterval()
    {
        return displayInterval;
    }


    public float getWellDirectionAngle()
    {
        if(targetPoint == null) return bigFloat;
        float deltaX = targetPoint.getX() - kickoffPoint.getX();
        float deltaY = targetPoint.getY() - kickoffPoint.getY();
        return StsMath.atan2(deltaY, deltaX);
    }

    /** The sequence of wellLegs begins with a startStraightLeg and ends with an endStraightLeg or inclinedStraightLeg.
     *  The wellLegs alternate between straightLegs and arcLegs.  A startArcLeg follows the startStraightLeg and a
     *  endArcLeg precedes the inclinedStraightLeg.
     *
     *  A build & hold consists of: startStraightLeg, startArcLeg, endStraightLeg.
     *  A build, hold, & drop consists of: startStraightLeg, startArcLeg, holdStraightLeg, endArcLeg, endStraightLeg.
     *
     *  If we insert an additional "dodge" point in a holdStraightLeg, we will have,
     *  if added to a build & hold: startStraightLeg, startArcLeg, holdLegOne, middleArcLeg, holdLegTwo, endStraightLeg;
     *  if added to a build, hold, & drop: startStraightLeg, startArcLeg, holdLegOne, middleArcLeg, holdLegTwo, endArcLeg, endStraightLeg.
     *
     *  Each straight leg has a vector from startPoint to endPoint; each arcLeg has an endPoint.  These vectors and points
     *  are referenced by connected legs so that each leg has an inPoint, outPoint, inVector, and outVector.
     *
     *  Certain legs have targetPoints which must be honored.  They are also used to setup an initial series of straight-line
     *  connections defining a rough approximation to the final path.  A startLeg has a targetPoint which is the kickoffPoint.  An
     *  inclinedLeg has a targetPoint which is the actual target.  When an itermediate "dodge" point is specified, this is a
     *  targetPoint for an arc which passes thru this point half-way thru the arc.
     *
     *  An inclinedTargetLeg is a straightLeg which descends at a given inclination from a height above the target and hits the target.
     *  This inclinedTargetLeg is proceeded by the endArcAngle which is in the vertical plane of the inclinedTargetLeg.  The azimuth
     *  of the inclinedTargetLeg generally points towards the last nonNull targetPoint which is a dodge point or the kickoffPoint.
     *
     *  After construction and setup, we iterate thru the sequence of legs until we have convergence.  Convergence is established when
     *  the center-point for all arcs are not moved from their position on the previous iteration.
     *
     *  Examine the compute() method for each subClass to see what is being computed during an iteration.
     *
     *  The test routine run from main() exercises a number of wellLeg construction scenarios.
     */

    abstract public class Leg
    {
        public Leg prevLeg = null;
        public Leg nextLeg = null;
        public String legName;
        public StsPoint inPoint;
        public StsPoint outPoint;
        public StsPoint inVector;
        public StsPoint outVector;
        public StsPoint targetPoint = null;
        public boolean converged = true;

        /** Automatically add each leg created to the legList.  If a new leg is to be inserted in between,
         *  it must be removed from the end of the list and inserted where required.
         */

        public Leg()
        {
            addLeg(this);
        }

        /** Search backwards for a targetPoint which can be used to shoot at. */
        public StsPoint getPrevTargetPoint()
        {
            if (targetPoint != null)
            {
                return targetPoint;
            }
            else
            {
                return prevLeg.getPrevTargetPoint();
            }
        }

        /** Search forwards for a targetPoint which can be used to shoot at. */
        public StsPoint getNextTargetPoint()
        {
            if (targetPoint != null)
            {
                return targetPoint;
            }
            else
            {
                return nextLeg.getNextTargetPoint();
            }
        }

        /** A prevLeg has asked this leg for its outPoint; if not available return
         *  this targetPoint or the next one that can be found.
         */

        public StsPoint getNextOutPoint()
        {
            if (outPoint.v[0] != bigFloat)
            {
                return outPoint;
            }
            else
            {
                return getNextTargetPoint();
            }
        }

        /** Used by an arcLeg wanting to define a vector from the prevInPoint of a straight leg to this inPoint. */
        public StsPoint getPrevInPoint()
        {
            if (!isPointNull(inPoint))
            {
                return inPoint;
            }
            else
            {
                return prevLeg.getPrevTargetPoint();
            }
        }

        /** A point has been created, but it hasn't been initialized yet: return true. */
        public boolean isPointNull(StsPoint point)
        {
            return point.v[0] == bigFloat;
        }

        /** Compute the geometry for this leg; override in subclass to actuall perform. */
        public void compute() // computes new leg
        {
        }

        /** adjust this leg iteratively; override in subclass as needed */
        public boolean adjust() // iterative adjustment of leg using neighbor legs
        {
            return false;
        }

        abstract public void initialize(); // initializes point and vector values

        abstract public void linkLegs(); // links point and vector values between legs

        abstract public void complete(); // computes time values for depth points

        abstract public void debug(); // mainDebug print out

        abstract public void display(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet);

        abstract public void pick(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet);

        abstract public StsPoint[] getPoints(); // gets leg points for assembly into planWell points

	    abstract public StsPoint[] getExportPoints(); // gets sufficient points for accurately exporting geometry
    }

    /** When an arc is inserted in an existing straight leg, insert the arc before the existing leg
     *  and then insert this new straight leg before the arc.
     */

    public StraightLeg constructStraightLegInsertBefore(String name, Leg nextLeg)
    {
        StsPoint outPoint = createNullPoint();
        MidStraightLeg wellLeg = new MidStraightLeg(name);
        wellLegsListInsertBefore(wellLeg, nextLeg);
        return wellLeg;
    }

    /** Create a point and assign first value to largFloat to indicate its not yet initialized. */
    private StsPoint createNullPoint()
    {
        StsPoint point = new StsPoint(5);
        point.v[0] = bigFloat;
        return point;
    }

    /** Add this wellLeg to the end of the current list of well legs. */
//    private void addToWellLegs(Leg wellLeg)
//    {
//        wellLegsList = (Leg[]) StsMath.arrayAddElement(wellLegsList, wellLeg);
//    }

    /** Insert this wellLeg before nextLeg in the list. */
    private void wellLegsListInsertBefore(Leg wellLeg, Leg nextLeg)
    {
        int nLegs = wellLegsList.length;
        int n;

        Leg matchLeg = null;
        for (n = 0; n < nLegs; n++)
        {
            matchLeg = wellLegsList[n];
            if (matchLeg == nextLeg)
            {
                break;
            }
        }
        if (matchLeg == null)
        {
            StsException.systemError("StsWellPlanSet.wellLegsInsertBefore() failed: could not find matching leg for " + nextLeg.legName);
            return;
        }
        // wellLeg has been automatically added to end of list, so delete it first and then insert it in proper location
        deleteLeg(wellLeg);
        insertLegAtIndex(wellLeg, n);
    }

    private void insertLegAtIndex(Leg wellLeg, int index)
    {
        wellLegsList = (Leg[]) StsMath.arrayInsertElementBefore(wellLegsList, wellLeg, index);
    }

    /** Abstract class for the 4 types of straightLegs: start, middle, end, and inclined. */
    abstract public class StraightLeg extends Leg
    {
        public StsPoint[] getPoints()
        {
            return new StsPoint[]
                {inPoint, outPoint};
        }

	    public StsPoint[] getExportPoints()
		{
			int nPoints = 8;
			int nIncrements = nPoints-1;
			StsPoint[] points = new StsPoint[nPoints];
			points[0] = new StsPoint(inPoint);
			for(int n = 1; n < nIncrements; n++)
				points[n] = StsPoint.staticInterpolatePoints(inPoint, outPoint, (float)n/nIncrements);
			points[nIncrements] = new StsPoint(outPoint);
			return points;
		}

        public void debug()
        {
            System.out.println(legName);
            System.out.print("    inPoint: ");
            inPoint.print();
            System.out.print("    outPoint: ");
            outPoint.print();
            if (targetPoint != null)
            {
                System.out.print("    targetPoint: ");
                targetPoint.print();
            }

            StsPoint vector = StsPoint.subPointsStatic(outPoint, inPoint);
            vector.normalize();
            System.out.print("    vector: ");
            vector.print();
        }

        public void complete()
        {
            computeTime(outPoint);
        }

        public void display(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet)
        {
            try
            {
                gl.glDisable(GL.GL_LIGHTING);
                if (highlighted)
                {
                    gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
                }
                else
                {
                    gl.glLineWidth(StsGraphicParameters.well3dLineWidth);
                }
                gl.glBegin(GL.GL_LINE_STRIP);

                color.setGLColor(gl);
                gl.glVertex3fv(inPoint.getPointXYZ(), 0);
                gl.glVertex3fv(outPoint.getPointXYZ(), 0);
            }
            catch (Exception e)
            {
                StsException.outputException("StsWellPlanStraightLeg.display() failed. ",
                                             e, StsException.WARNING);
            }
            finally
            {
                gl.glEnd();
                gl.glEnable(GL.GL_LIGHTING);
            }
        }

        public void pick(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet)
        {
            try
            {
                gl.glDisable(GL.GL_LIGHTING);
                if (highlighted)
                {
                    gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
                }
                else
                {
                    gl.glLineWidth(StsGraphicParameters.well3dLineWidth);
                }
                gl.glBegin(GL.GL_LINE_STRIP);

                color.setGLColor(gl);
                gl.glVertex3fv(inPoint.getPointXYZ(), 0);
                gl.glVertex3fv(outPoint.getPointXYZ(), 0);
            }
            catch (Exception e)
            {
                StsException.outputException("StsWellPlanStraightLeg.display() failed. ",
                                             e, StsException.WARNING);
            }
            finally
            {
                gl.glEnd();
                gl.glEnable(GL.GL_LIGHTING);
            }
        }
    }

    /** The first down leg from kb to kickoff. */
    public class StartStraightLeg extends StraightLeg
    {
        public StartStraightLeg(String name, StsPoint kbPoint, StsPoint kickoffPoint)
        {
            this.legName = name;
            inPoint = kbPoint;
            outPoint = kickoffPoint;
        }

        public void initialize()
        {
            inVector = StsPoint.subPointsStatic(outPoint, inPoint);
            inVector.normalize();
            outVector = inVector;
            targetPoint = outPoint;
        }

        public void linkLegs()
        {
        }
    }

    /** a middleLeg between two arcs. */
    public class MidStraightLeg extends StraightLeg
    {
        public MidStraightLeg(String name)
        {
            legName = name;
        }

        public void initialize()
        {
            outVector = createNullPoint();
            inVector = outVector;
            outPoint = createNullPoint();
        }

        public void linkLegs()
        {
            inPoint = prevLeg.outPoint;
            prevLeg.outVector = inVector;
        }
    }

    /** An arbitary endLeg which has not fixed inclination. */
    public class EndStraightLeg extends StraightLeg
    {
        public EndStraightLeg(String name, StsPoint targetPoint)
        {
            legName = name;
            this.targetPoint = new StsPoint(targetPoint);
        }

        public void initialize()
        {
            if(outVector == null) outVector = createNullPoint();
            if(inVector == null) inVector = outVector;
            outPoint = targetPoint;
        }
/*

        public void initializeInVector(StsPoint inPoint)
        {
           inVector = StsPoint.subPointsStatic(targetPoint, inPoint);
           inVector.normalize();
        }
*/

        public void linkLegs()
        {
            inPoint = prevLeg.outPoint;
            prevLeg.outVector = inVector;
        }
    }

    /** An inclined leg specified by inclination of vertical and vertical distance from top of straight leg
     *  down to target.
     */

    public class InclinedEndLeg extends EndStraightLeg
    {
        public float dropoffInclinationAngle;
        public float dzAtDropoffEndAboveTarget;

        public InclinedEndLeg(String name, float dropoffInclinationAngle, float dzAtDropoffEndAboveTarget, StsPoint targetPoint)
        {
            super(name, targetPoint);
            this.dropoffInclinationAngle = dropoffInclinationAngle;
            this.dzAtDropoffEndAboveTarget = dzAtDropoffEndAboveTarget;
        }

        public void compute()
        {
            float dx, dy;

            StsPoint prevInPoint = prevLeg.getPrevInPoint();
            dx = targetPoint.getX() - prevInPoint.getX();
            dy = targetPoint.getY() - prevInPoint.getY();

            StsPoint horizontalVector = new StsPoint(5);
            horizontalVector.setX(dx);
            horizontalVector.setY(dy);
            horizontalVector.setZ(0.0f);
            horizontalVector.normalize();
            verticalVector.setX(0.0f);
            verticalVector.setY(0.0f);
            verticalVector.setZ(1.0f);
            float tanAngle = (float) Math.tan(dropoffInclinationAngle * StsMath.RADperDEG);
            inVector.multByConstantAddPoint(horizontalVector, tanAngle, verticalVector);
            // Here vertical component of inVector is length 1.0 and horizontal component is length is tanAngle.
            // We want vertical component to have length dzAtDropoffEndAboveTarget, so scale by this times -1
            // and add to outPoint (targetPoint) to get the inPoint.
            inPoint.multByConstantAddPoint(inVector, -dzAtDropoffEndAboveTarget, outPoint);
            inVector.normalize();
        }

        public boolean adjust()
        {
            StsPoint oldInPoint = new StsPoint(inPoint);
            compute();
            converged = inPoint.sameAs(oldInPoint);
            return !converged;
        }

        public void debug()
        {
            System.out.println(legName);
            System.out.print("    inPoint: ");
            inPoint.print();
            System.out.print("    outPoint: ");
            outPoint.print();
            float dx = inVector.getX();
            float dy = inVector.getY();
            double azimuth = StsMath.DEGperRAD * Math.atan2(dy, dx);
            System.out.println("    azimuth: " + azimuth);
        }
    }

    /** Abstract class for the 3 types of arcLegs: start, middle, and end. */
    abstract public class ArcLeg extends Leg
    {
        public float radius;
        public float adjustedRadius;
        StsPoint axis;
        float angleRad;
        StsPoint center;
        StsPoint[] points = null;
        int nIncrements = 0;
        float maxError = 0.0001f;

        static final long serialVersionUID = 1L;

        public ArcLeg(String name, float radius)
        {
            super();
            this.legName = name;
//            this.type = type;
            this.radius = radius;
//            this.adjustedRadius = radius;
//            outPoint = createNullPoint();
//            outVector = createNullPoint();
        }

        abstract public void initialize();

        /** prevLeg must be a straightLeg.  Construct right triangle from center to prevInPoint
         *  and back to tangent point which we want to find.  If center to inPoint is less than
         *  radius of arc, we need to reduce the radius to something just less than this so we
         *  can find a solution.  This implies radius is too tight.  Once we have a solution,
         *  may need to check if this constraint can be slowly loosened up.
         */
        public void computeInPoint()
        {
            float rotAngleRad;
            StsPoint prevInPoint = prevLeg.inPoint;
            StsPoint centerToPrevInPointVector = StsPoint.subPointsStatic(prevInPoint, center);
            float length = centerToPrevInPointVector.normalizeXYZReturnLength();
            if (length < adjustedRadius)
            {
                System.out.print("radius for " + this.legName + " adjusted down from " + adjustedRadius + " to ");
                adjustedRadius = 0.9f * length;
                System.out.println(adjustedRadius);
            }
            /*
                else if(adjustedRadius < radius && length > adjustedRadius/0.98f)
                {
                    System.out.print("radius for " + this.legName + " adjusted up from " + adjustedRadius + " to ");
                    float dif = Math.min(radius, length) - adjustedRadius;
                    adjustedRadius /= 0.99f;
                    System.out.println(adjustedRadius);
                }
             */
            rotAngleRad = (float) Math.acos(adjustedRadius / length);
            centerToPrevInPointVector.multiply(adjustedRadius);
            rotatePointAroundAxis(centerToPrevInPointVector, center, axis, rotAngleRad, inPoint);
            inVector.subPoints(inPoint, prevInPoint);
            inVector.normalize();
            // inPoint and inVector have changed, so recompute
            axis = StsPoint.leftCrossProductStatic(inVector, outVector);
            length = (float) axis.length();
            float dot = inVector.dot(outVector, 3);
            angleRad = (float) Math.atan2(length, dot);
        }

        /** Point is defined as a vector and origin. Rotation is around axis thru angleRad. Result is stored in rotatedPoint. */
        private void rotatePointAroundAxis(StsPoint vector, StsPoint origin, StsPoint axis, float angleRad, StsPoint rotatedPoint)
        {
            StsQnion q = new StsQnion(axis, angleRad);
            StsPoint rotatedVector = q.leftRotateVec(vector);
            rotatedPoint.addPoints(origin, rotatedVector);
        }

        /** Same geometric algorithm as computeInPoint() but applied to the outVector. */
        public void computeOutPoint()
        {
            float rotAngleRad;
            StsPoint nextOutPoint = nextLeg.getNextOutPoint();
            StsPoint centerToNextOutPointVector = StsPoint.subPointsStatic(nextOutPoint, center);
            float length = centerToNextOutPointVector.normalizeXYZReturnLength();
            if (length < adjustedRadius)
            {
                System.out.print("radius for " + this.legName + " adjusted down from " + adjustedRadius + " to ");
                adjustedRadius = 0.9f * length;
                System.out.println(adjustedRadius);
            }
            /*
                else if(adjustedRadius < radius && length > adjustedRadius/0.98f)
                {
                    System.out.print("radius for " + this.legName + " adjusted up from " + adjustedRadius + " to ");
                    float dif = Math.min(radius, length) - adjustedRadius;
                    adjustedRadius /= 0.99f;
                    System.out.println(adjustedRadius);
                }
             */
            rotAngleRad = (float) Math.acos(adjustedRadius / length);
            centerToNextOutPointVector.multiply(adjustedRadius);
            rotatePointAroundAxis(centerToNextOutPointVector, center, axis, -rotAngleRad, outPoint);
            outVector.subPoints(nextOutPoint, outPoint);
            outVector.normalize();
            axis = StsPoint.leftCrossProductStatic(inVector, outVector);
            length = axis.length();
            float dot = inVector.dot(outVector, 3);
            angleRad = (float) Math.atan2(length, dot);
        }

        public void approximateOutVector()
        {
            StsPoint nextOutPoint = nextLeg.getNextOutPoint();
            outVector.subPoints(nextOutPoint, inPoint);
            outVector.normalize();
        }

        public void computeCenterOffInPoint()
        {
            if(isPointNull(outVector)) approximateOutVector();
            axis = StsPoint.leftCrossProductStatic(inVector, outVector);
            float length = axis.length();
            float dot = inVector.dot(outVector, 3);
            angleRad = (float) Math.atan2(length, dot);
            StsPoint inRadialVector = StsPoint.leftCrossProductStatic(inVector, axis); // center to inPoint vector
            inRadialVector.normalize();
            inRadialVector.multiply(adjustedRadius);
            center = StsPoint.subPointsStatic(inPoint, inRadialVector);
        }

        public boolean adjust()
        {
            converged = true;
            StsPoint oldCenter = null;
            if (center != null)
            {
                oldCenter = new StsPoint(center);
            }
            StsPoint oldAxis = null;
            if (axis != null)
            {
                oldAxis = new StsPoint(axis);
            }
            compute();
            if (oldCenter == null && oldAxis == null)
            {
                converged = false;
            }
            if (!oldCenter.sameAs(center))
            {
                converged = false;
            }
            if (!oldAxis.sameAs(axis))
            {
                converged = false;
            }
            return !converged;
        }

        public StsPoint[] getPoints()
        {
            computePoints();
            return points;
        }

	    public StsPoint[] getExportPoints() { return getPoints(); }

        public void complete()
        {
//            computeTime(inPoint);
            computeTime(outPoint);
//            computeCenterPoint();
            if (debug)
            {
                if (targetPoint != null)
                {
                    computeTime(targetPoint);
                }
                if (adjustedRadius < radius)
                {
                    System.out.println("RADIUS reduced for " + legName + " from " + radius + " to " + adjustedRadius);
                }
            }
        }

        public void display(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet)
        {
            if (points == null) computePoints();
            StsGLDraw.drawLine(gl, color, highlighted, points, zIndex);

        }

        public void pick(GL gl, StsColor color, boolean highlighted, int zIndex, StsWellPlan plannedWellSet)
        {
            if (points == null) computePoints();
            StsGLDraw.pickLine(gl, color, highlighted, points, zIndex);
        }

        private StsPoint[] computePoints()
        {
            int nIncrements = Math.max(10, StsMath.above(angleRad * adjustedRadius / displayInterval));
            double dAngle = angleRad / nIncrements;
            int nPoints = nIncrements + 1;
            points = new StsPoint[nPoints];
            points[0] = inPoint;
            StsPoint radial = StsPoint.subPointsStatic(inPoint, center);
            StsPoint rotator = new StsPoint(axis);
            rotator.normalize();
            float tanAngle = (float) Math.tan(dAngle);
            rotator.multiply(tanAngle);
            float cosAngle = (float) Math.cos(dAngle);
            for (int n = 1; n < nPoints; n++)
            {
                StsPoint dTangent = StsPoint.leftCrossProductStatic(rotator, radial);
                radial.add(dTangent);
                radial.multiply(cosAngle);
                StsPoint radialPoint = StsPoint.addPointsStatic(radial, center);
                computeTime(radialPoint);
                points[n] = radialPoint;
            }
            //           points[nPoints - 1] = outPoint;
            return points;
        }

        public void debug()
        {
            System.out.println(legName);
            System.out.print("    inPoint: ");
            inPoint.print();
            System.out.println("    radius: " + adjustedRadius);
            float turnRate = (float) (18000 / (Math.PI * adjustedRadius));
            System.out.println("    turn rate: " + turnRate + " deg/100");
            System.out.println("    angle: " + angleRad * StsMath.DEGperRAD);
            System.out.print("    inVector: ");
            inVector.print();
            if (axis != null)
            {
                System.out.print("    axis: ");
                axis.print();
            }
            if (center != null)
            {
                System.out.print("    center: ");
                center.print();
            }
            System.out.print("    outPoint: ");
            outPoint.print();
            System.out.print("    outVector: ");
            outVector.print();
            if (targetPoint != null)
            {
                System.out.print("    targetPoint: ");
                targetPoint.print();
            }
            if (axis != null)
            {
                StsPoint[] points = computePoints();
                int nPoints = points.length;
                System.out.println("Arc points");
                for (int n = 0; n < nPoints; n++)
                {
                    System.out.print("   n " + n + " ");
                    points[n].print();
                }
            }
        }
    }

    /** The arc immediately after the StartStraightLeg.  Begins at the kickoffPoint
     *  and proceeds along arc in direction of next targetPoint.
     */
    public class StartArcLeg extends ArcLeg
    {
        public StartArcLeg(String name, float radius)
        {
            super(name, radius);
        }

        public void initialize()
        {
            outPoint = createNullPoint();
            adjustedRadius = radius;
        }

        public void linkLegs()
        {
            inPoint = prevLeg.outPoint;
            inVector = prevLeg.outVector;
            outVector = nextLeg.inVector;
        }

        public void compute()
        {
            computeCenterOffInPoint();
            computeOutPoint();
        }
    }

    /** The midArcLeg is between two midStraightLegs. */
    public class MidArcLeg extends ArcLeg
    {
        public MidArcLeg(String name, float radius, StsPoint midTargetPoint)
        {
            super(name, radius);
            targetPoint = new StsPoint(midTargetPoint);
//            outPoint = new StsPoint(midTargetPoint);
        }

        public void initialize()
        {
            outPoint = new StsPoint(targetPoint);
            adjustedRadius = radius;
        }

        public void linkLegs()
        {
            inVector = prevLeg.outVector;
            outVector = nextLeg.inVector;
            inPoint = prevLeg.outPoint;
            inPoint.copyFrom(targetPoint);
        }

        public void compute()
        {
            computeCenterPoint();
            computeInPoint();
            computeOutPoint();
        }

        public boolean adjust()
        {
            converged = true;
            StsPoint oldCenter = null;
            if (center != null)
            {
                oldCenter = new StsPoint(center);
            }
            compute();
            if (oldCenter == null)
            {
                converged = false;
            }
            converged = oldCenter.sameAs(center);
            if (debug)
            {
                if (!converged)
                {
                    System.out.println("    for " + legName + " center moved from " + oldCenter.toString() + " to " + center.toString());
                }
                else
                {
                    System.out.println("    midArc " + legName + " converged.");
                }
            }
            return !converged;
        }

        public void computeCenterPoint()
        {
            axis = StsPoint.leftCrossProductStatic(inVector, outVector);
            StsPoint centerVector = StsPoint.subPointsStatic(outVector, inVector);
            centerVector.normalize();
            center = StsPoint.multByConstantAddPointStatic(centerVector, adjustedRadius, targetPoint);
        }
    }

    /** The arcLeg just before the inclinedEndLeg. */
    public class EndArcLeg extends ArcLeg
    {
        public EndArcLeg(String name, float radius)
        {
            super(name, radius);
        }

        public void initialize()
        {
            outPoint = createNullPoint();
            adjustedRadius = radius;
        }

        public void linkLegs()
        {
            inVector = prevLeg.outVector;
            outVector = nextLeg.inVector;
            inPoint = prevLeg.outPoint;
        }

        public void compute()
        {
            // compute location of center
            axis = StsPoint.leftCrossProductStatic(outVector, verticalVector);
            if (axis.length() < StsMath.FLT_EPSILON)
            {
                axis = StsPoint.leftCrossProductStatic(inVector, verticalVector);
            }
            StsPoint outRadialVector = StsPoint.leftCrossProductStatic(axis, outVector);
            outRadialVector.normalize();
            outRadialVector.multiply(adjustedRadius);
            center = StsPoint.addPointsStatic(outPoint, outRadialVector);
            computeInPoint();
        }
    }

	/** The first down leg from kb to kickoff. */
	 public class StartWellArcLeg extends ArcLeg
	 {
		 public StartWellArcLeg(String name, float radius, StsPoint inPoint, StsPoint inVector)
		 {
			 super(name, radius);
			 this.inPoint = new StsPoint(inPoint);
			 this.inVector = inVector;
		 }

		 public void initialize()
		 {
			 outPoint = createNullPoint();
			 adjustedRadius = radius;
		 }

         public void linkLegs()
		 {
			 outVector = nextLeg.inVector;
		 }

		 public void compute()
		 {
			 computeCenterOffInPoint();
			 computeOutPoint();
		 }
    }

    public void debugTest()
    {
        StraightLeg kbLeg;
        InclinedEndLeg targetLeg;
        ArcLeg kickoffArcLeg;
        StraightLeg holdLeg;
        ArcLeg dodgeArcLeg1, dodgeArcLeg2;
        StraightLeg holdLeg1, holdLeg2, holdLeg3;
        ArcLeg dropArcLeg;
		StartWellArcLeg startWellLeg;

        StsPoint kbPoint = new StsPoint(0.0f, 100.0f, -10.0f);
        StsPoint kickoffPoint = new StsPoint(0.0f, 100.0f, 0.0f);
        StsPoint targetPoint = new StsPoint(50f, 150.0f, 50.0f);
        float buildRadius = 10.0f;
        float midRadius = 10.0f;

        System.out.println("\nBUILD AND HOLD");
        initializeWellLegsList();
        kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
        kickoffArcLeg = new StartArcLeg("Kickoff arc", buildRadius);
        holdLeg = new EndStraightLeg("Hold leg", targetPoint);
        doLegs();

        System.out.println("\nBUILD, HOLD, AND ADD DODGE");
        StsPoint dodgePoint = new StsPoint(12.0f, 100.0f, 12.0f);
        insertArcLegInStraightLeg("Dodge arc 1", midRadius, dodgePoint, holdLeg);
        doLegs();

//        targetPoint = new StsPoint(30.f, 100.0f, 40.0f);
        float dropoffInclinationAngle = 30.0f;
        float dzAtDropoffEndAboveTarget = 10.0f;
        float dropRadius = 10.0f;
        StsPoint dodgePoint1 = new StsPoint(30.0f, 100.0f, 0.0f);
        StsPoint dodgePoint2 = new StsPoint(20.0f, 100.0f, 40.0f);

        System.out.println("\nBUILD, HOLD AND DROP");

        initializeWellLegsList();
        kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
        kickoffArcLeg = new StartArcLeg("Kickoff arc", buildRadius);
        holdLeg = new MidStraightLeg("Hold leg");
        dropArcLeg = new EndArcLeg("Drop arc", dropRadius);
        targetLeg = new InclinedEndLeg("Inclined target leg", dropoffInclinationAngle, dzAtDropoffEndAboveTarget, targetPoint);
        doLegs();

        System.out.println("\nBUILD, HOLD, DROP AND ADD DODGE");
        insertArcLegInStraightLeg("Dodge arc 1", midRadius, dodgePoint1, holdLeg);
        doLegs();

        System.out.println("\nBUILD, HOLD, DROP AND ADD TWO DODGES");
        insertArcLegInStraightLeg("Dodge arc 2", midRadius, dodgePoint2, holdLeg);
        doLegs();

        System.out.println("\nBUILD, HOLD TWO DODGES AND DROP");

        initializeWellLegsList();
        kbLeg = new StartStraightLeg("KB leg", kbPoint, kickoffPoint); // defines kb to ko leg
        kickoffArcLeg = new StartArcLeg("Kickoff arc", buildRadius);
        holdLeg1 = new MidStraightLeg("Hold leg 1");
        dodgeArcLeg1 = new MidArcLeg("Dodge arc 1", midRadius, dodgePoint1);
        holdLeg2 = new MidStraightLeg("Hold leg 2");
        dodgeArcLeg2 = new MidArcLeg("Dodge arc 2", midRadius, dodgePoint2);
        holdLeg3 = new MidStraightLeg("Hold leg 3");
        dropArcLeg = new EndArcLeg("Drop arc", dropRadius);
        targetLeg = new InclinedEndLeg("Inclined target leg", dropoffInclinationAngle, dzAtDropoffEndAboveTarget, targetPoint);
        /*
                 kickoffArcLeg.radius = 7.3f;
                 dodgeArcLeg1.radius = 8.1f;
                 dodgeArcLeg2.radius = 1.2f;
                 dropArcLeg.radius = 6.8f;
                 kickoffArcLeg.adjustedRadius = 7.3f;
                 dodgeArcLeg1.adjustedRadius = 8.1f;
                 dodgeArcLeg2.adjustedRadius = 1.2f;
                 dropArcLeg.adjustedRadius = 6.8f;
         */
        doLegs();

        System.out.println("\nBUILD, HOLD TWO DODGES AND DROP, DELETE FIRST DODGE");
        deleteLegs(new Leg[]
                   {holdLeg2, dodgeArcLeg2});
        doLegs();

		System.out.println("\nEXTEND WELL WITH HOLD LEG");
	    initializeWellLegsList();
		StsPoint startPoint = new StsPoint(0.0f, 0.0f, 0.0f);
		StsPoint startVector = new StsPoint(1.0f, 0.0f, 1.0f);
		startWellLeg = new StartWellArcLeg("Start well Leg", buildRadius, startPoint, startVector);
		holdLeg = new EndStraightLeg("Hold leg", targetPoint);
		doLegs();
    }


	public StsSurfaceVertex addAbsolutePoint(double x, double y, double z, double m, double t)
	{
		StsPoint point = new StsPoint(5);
		point.setX((float)(x - xOrigin));
		point.setY((float)(y - yOrigin));
		point.setZ((float)z);
		point.setT((float)t);
		point.setM((float)m);
		StsSurfaceVertex vertex = addLineVertex(point, true, false);
		computePoints();
		currentModel.win3dDisplayAll();
		return vertex;
	}

	public StsSurfaceVertex addAbsolutePoint(double x, double y, double z)
	{
		StsPoint point = new StsPoint(5);
		point.setX((float)(x - xOrigin));
		point.setY((float)(y - yOrigin));
		point.setZ((float)z);
		StsSurfaceVertex vertex = addLineVertex(point, true, false);
		computePoints();
		currentModel.win3dDisplayAll();
		return vertex;
	}

	public void setAbsoluteXYZ(double x, double y, double z)
	{
		StsPoint point = currentVertex.getPoint();
		point.setX((float)(x - xOrigin));
		point.setY((float)(y - yOrigin));
		point.setZorT((float)z);
		computePoints();
		currentModel.win3dDisplayAll();
	}

	public double getAbsoluteX()
	{
		if(currentVertex == null)
			return 0.0f;
		double[] xy = project.getAbsoluteXYCoordinates(currentVertex.getPoint());
		return xy[0];
	}

	public double getAbsoluteY()
	{
		if(currentVertex == null)
			return 0.0f;
		double[] xy = project.getAbsoluteXYCoordinates(currentVertex.getPoint());
		return xy[1];
	}

	public double getAbsoluteZ()
	{
		if(currentVertex == null)
			return 0.0f;
		return currentVertex.getPoint().getZorT();
	}

	public StsSurfaceVertex addLineVertex(StsPoint point)
	{
		StsSurfaceVertex vertex = super.addLineVertex(point, false, false);
		if(vertex == null) return null;

		setCurrentVertex(vertex);
		computePoints();
		return vertex;
	}

	public int getVertexIndex(StsSurfaceVertex vertex)
	{
		if(lineVertices == null) return -1;
		return lineVertices.getIndex(vertex);
	}

	public void checkChangeCurrentVertex()
	{
		setCurrentVertex(getTopLineVertex());
	}

	public void setCurrentVertex(StsSurfaceVertex vertex)
	{
		currentVertex = vertex;
		setSelectedVertex(vertex);
	}
	public StsSurfaceVertex getCurrentVertex() { return currentVertex; }

	public void pickVertices(StsGLPanel3d glPanel3d)
	{
		super.pickVertices(glPanel3d);
	}


	public boolean buildWell(StsWell actualWell)
	{
		// Add the actual well points to the list.
		if(actualWell == null)
			return false;
		setXOrigin(actualWell.xOrigin);
		setYOrigin(actualWell.yOrigin);
		StsSurfaceVertex lastVertex = (StsSurfaceVertex)actualWell.getLineVertices().getLast();
		lineVertices.add(lastVertex);
		setCurrentVertex(lastVertex);
		computePoints();
		return true;
    }

	public boolean exportPlan(String name)
	{
		 try
		 {
			 String directory = currentModel.getProject().getDataFullDirString();
			 String filename = "wellPlan.obj." + name;
			 StsDBFileObjectTrader.exportStsObject(directory + File.separator + filename, this, null);

             //StsDBFile.writeObjectFile(directory + File.separator + filename, this, null);
//             StsToolkit.serializeObject(this, directory, filename);
			 return true;
		 }
		 catch(Exception e)
		 {
			 StsException.outputException("StsWellPlanExport() failed for " + name, e, StsException.WARNING);
			 return false;
		 }
	}

    public boolean exportPlanWell(String timeOrDepth)
    {
        return StsPlanWellExportDialog.exportWell(currentModel, currentModel.win3d, "Well Export Utility", true, this, timeOrDepth);
    }


    public boolean exportPath(String timeOrDepth)
	{
		if(drillingWell == null) return false;
		return drillingWell.export(timeOrDepth);
	}

    public static void main(String[] args)
    {

        StsWellPlan plannedWellSet = new StsWellPlan("test", false);
		plannedWellSet.initializePlanPoints();
        plannedWellSet.debug = true;
        plannedWellSet.debugTest();
    }

    public Class[] getViewClasses()
    {
        return new Class[] { StsSeismicCurtainView.class };
    }
}