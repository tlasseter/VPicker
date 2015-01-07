package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.UI.Beans.*;

//public class StsSeismicVelocityModelClass extends StsClass implements DBSerializable
public class StsPreStackVelocityModelClass extends StsClass implements StsSerializable
{
    protected boolean displayAnalysisPoints = true;

    public StsPreStackVelocityModelClass()
    {
        userName = "Pre-Stack Velocity Model";
    }

    public void initializeDisplayFields()
    {
        StsBooleanFieldBean displayAnalysisPointsBean = new StsBooleanFieldBean(this, "displayAnalysisPoints", "Display analysis points");
        displayFields = new StsFieldBean[]{displayAnalysisPointsBean};
    }

    public void setDisplayAnalysisPoints(boolean display) { displayAnalysisPoints = display; }

    public boolean getDisplayAnalysisPoints() { return displayAnalysisPoints; }
}