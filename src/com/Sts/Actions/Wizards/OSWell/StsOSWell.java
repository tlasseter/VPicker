package com.Sts.Actions.Wizards.OSWell;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;
import com.openspirit.data.*;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 10, 2008
 * Time: 7:53:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsOSWell
{
	// basic initial attributes for selection
    private String name = null;
    private String id = null;
    private double xOrigin = 0.0;
    private double yOrigin = 0.0;
    private DataKey dataKey = null;

    // data attributes for loading
    private String operator = "unknown";
    private String company = "unknown";
    private String field = "unknown";
    private String wellNumber = "unknown";
    private String wellLabel = "unknown";
    private String api = "00000000000000";
    private String uwi = "00000000000000";
    private String date = "unknown";
    private float kbElev = 0.f;
    private float elev = 0.f;
    private String elevDatum = "unknown";
    private long completionDate = 0L;
    private long spudDate = 0L;

    // deviation data
    private StsFloatVector devMD = null;
    private StsFloatVector devTVD = null;
    private StsFloatVector devTWT = null;
    private StsFloatVector devXOffset = null;
    private StsFloatVector devYOffset = null;

    private ArrayList<StsOSWellMarker> markers = null;
    private ArrayList<StsLogCurve> logs = null;

    private StsLogCurve tdCurve = null;

    /**
     * Initial constructor for list wells.
     * @param id is the well identifier
     * @param name is the common well name
     * @param xOrigin is the surface X coordinate
     * @param yOrigin is the surface Y coordinate
     * @param dataKey is the OpenSpirit PrimaryKey$ for this well
     */
    public StsOSWell(String id, String name, double xOrigin,
    		double yOrigin, DataKey dataKey)
    {
    	this.id = id;
        this.name = name;
        this.xOrigin= xOrigin;
        this.yOrigin = yOrigin;
        this.dataKey = dataKey;
    }

	public DataKey getDataKey()
	{
		return dataKey;
	}

	public void setDataKey(DataKey dataKey)
	{
		this.dataKey = dataKey;
	}

	public String getId()
	{
		return id;
	}

	public void setId(String id)
	{
		this.id = id;
	}

	public String getName()
	{
		return name;
	}

	public void setName(String name)
	{
		this.name = name;
	}

	public double getxOrigin()
	{
		return xOrigin;
	}

	public void setxOrigin(double xOrigin)
	{
		this.xOrigin = xOrigin;
	}

	public double getyOrigin()
	{
		return yOrigin;
	}

	public void setyOrigin(double yOrigin)
	{
		this.yOrigin = yOrigin;
	}

    public String toString() { return name; }

	public String getApi()
	{
		return api;
	}

	public void setApi(String api)
	{
		this.api = api;
	}

	public String getCompany()
	{
		return company;
	}

	public void setCompany(String company)
	{
		this.company = company;
	}

	public String getDate()
	{
		return date;
	}

	public void setDate(String date)
	{
		this.date = date;
	}

	public String getField()
	{
		return field;
	}

	public void setField(String field)
	{
		this.field = field;
	}

	public float getKbElev()
	{
		return kbElev;
	}

	public void setKbElev(float kbElev)
	{
		this.kbElev = kbElev;
	}

    public float getElev()
    {
        return elev;
    }

    public void setElev(float elev)
    {
        this.elev = elev;
    }

    public long getCompletionDate()
    {
        return completionDate;
    }

    public void setCompletionDate(long complete)
    {
        this.completionDate = complete;
    }

    public long getSpudDate()
    {
        return spudDate;
    }

    public void setSpudDate(long spud)
    {
        this.spudDate = spud;
    }

    public String getElevDatum()
    {
        return elevDatum;
    }

    public void setElevDatum(String datum)
    {
        this.elevDatum = datum;
	}
	public String getOperator()
	{
		return operator;
	}

	public void setOperator(String operator)
	{
		this.operator = operator;
	}

	public String getUwi()
	{
		return uwi;
	}

	public void setUwi(String uwi)
	{
		this.uwi = uwi;
	}

	public String getWellLabel()
	{
		return wellLabel;
	}

	public void setWellLabel(String wellLabel)
	{
		this.wellLabel = wellLabel;
	}

	public void setDevData(StsFloatVector md, StsFloatVector tvd, StsFloatVector twt,
			StsFloatVector xoffset, StsFloatVector yoffset)
	{
		devMD = md;
		devTVD = tvd;
		devTWT = twt;
		devXOffset = xoffset;
		devYOffset = yoffset;
	}

	public StsFloatVector getDevMD()
	{
		return devMD;
	}

	public StsFloatVector getDevTVD()
	{
		return devTVD;
	}

	public StsFloatVector getDevTWT()
	{
		return devTWT;
	}

	public StsFloatVector getDevXOffset()
	{
		return devXOffset;
	}

	public StsFloatVector getDevYOffset()
	{
		return devYOffset;
	}

	public void addWellMarker(StsOSWellMarker marker)
	{
		if (markers == null)
			markers = new ArrayList<StsOSWellMarker>();

		markers.add(marker);
	}

	public StsOSWellMarker[] getMarkers()
	{
		if (markers == null) return null;
		StsOSWellMarker[] mrkers = new StsOSWellMarker[markers.size()];
		markers.toArray(mrkers);
		return mrkers;
	}

	public void addLog(StsLogCurve log)
	{
		if (logs == null)
			logs = new ArrayList<StsLogCurve>();

		logs.add(log);
	}

	public StsLogCurve[] getLogs()
	{
		if (logs == null) return null;
		StsLogCurve[] curves = new StsLogCurve[logs.size()];
		logs.toArray(curves);
		return curves;
	}

	public void setTdCurve(StsLogCurve tdCurve)
	{
		this.tdCurve = tdCurve;
	}

	public StsLogCurve getTdCurve()
	{
		return tdCurve;
	}

}
