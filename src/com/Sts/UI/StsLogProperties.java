package com.Sts.UI;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;


public class StsLogProperties extends StsPanelProperties
{
	transient StsGroupBox groupBox = null;

    double datum = 0.0f;
    double velocity = 1000.0f;
    int dAttributeIndex = 0, vAttributeIndex = 0;
    boolean timeDatum = false;

    transient String[] dAttributes = null;
    transient String[] vAttributes = null;
    transient StsBooleanFieldBean timeBtn = new StsBooleanFieldBean();

    transient StsComboBoxFieldBean datumSelectBean = null;
    transient StsComboBoxFieldBean velocitySelectBean = null;
    transient StsDoubleFieldBean datumBean;
    transient StsDoubleFieldBean velocityBean;

    transient static public final byte NONE = 0;
    transient static public final byte USER_SPECIFIED = 1;
    transient static public String[] attributeStrings = new String[] {"None", "User Specified"};

	transient public StsModel model;

	transient public boolean recompute = true;

	static private final String title = "Log Properties";

	public StsLogProperties()
	{
	}

	public StsLogProperties(StsModel model, StsClass stsClass, String fieldName)
	{
		super(title, fieldName);
		this.model = model;
	}

	public StsLogProperties(StsObject parentObject, StsLogProperties defaultProperties, String fieldName)
	{
        super(parentObject, title, fieldName);
        initializeDefaultProperties(defaultProperties);
    }
    public void initializeDefaultProperties(Object defaultProperties)
    {
        if(defaultProperties != null)
            StsToolkit.copyObjectNonTransientFields(defaultProperties, this);
    }

	public void buildBeanLists(Object parentObject)
	{
        if(parentObject == null)
            return;
        String[] datumAttributes = null;
        String[] velocityAttributes = null;
        vAttributes = new String[] {"User Specified"};
        getDatumAttributes();

        if(parentObject instanceof StsPreStackLineSet)
        {
            velocityAttributes = ((StsPreStackLineSet)parentObject).lines[0].getVelocityAttributes();
            if(velocityAttributes != null)vAttributes = (String[])StsMath.arrayAddArray(vAttributes, velocityAttributes);
        }
        else if(parentObject instanceof StsVsp)
        {
            velocityAttributes = ((StsVsp)parentObject).getVelocityAttributes();
            if(velocityAttributes != null)vAttributes = (String[])StsMath.arrayAddArray(vAttributes, velocityAttributes);
		}
	}

    private void getDatumAttributes()
    {
        String[] datumAttributes = null;

        if(parentObject instanceof StsPreStackLineSet)
        {
           if(timeDatum)
                datumAttributes = ((StsPreStackLineSet) parentObject).lines[0].getTimeAttributes();
           else
                datumAttributes = ((StsPreStackLineSet) parentObject).lines[0].getDistanceAttributes();
        }
        else if(parentObject instanceof StsVsp)
        {
            if(timeDatum)
                datumAttributes = ((StsVsp) parentObject).getTimeAttributes();
            else
                datumAttributes = ((StsVsp) parentObject).getDistanceAttributes();
        }
        dAttributes = attributeStrings;
        if(datumAttributes != null)
            dAttributes = (String[])StsMath.arrayAddArray(dAttributes, datumAttributes);
    }

    public void initializeBeans()
    {
        buildBeanLists(parentObject);
        velocitySelectBean = new StsComboBoxFieldBean(this, "velocityAttribute", "Velocity for Datum Correction:", vAttributes);
        datumSelectBean = new StsComboBoxFieldBean(this, "datumAttribute", "Correct Traces to Datum:", dAttributes);
        timeBtn.initialize(this, "isTimeDatum", "Time Datum:");
        propertyBeans = new StsFieldBean[]
        {
           timeBtn,
           datumSelectBean,
           datumBean = new StsDoubleFieldBean(this, "datum", -10000, 10000, "Datum (m-ft):", false),
           velocitySelectBean,
           velocityBean = new StsDoubleFieldBean(this, "velocity", -10000, 10000, "Velocity (m-ft/sec):", false)
        };
        reconfigureUI();
    }

    public void setIsTimeDatum(boolean value)
    {
        timeDatum = value;
        getDatumAttributes();
        datumSelectBean.setListItems(dAttributes);
        if(dAttributeIndex > 1)
            dAttributeIndex = 0;
        datumSelectBean.setSelectedIndex(dAttributeIndex);
        velocitySelectBean.setSelectedIndex(vAttributeIndex);
        reconfigureUI();
    }
    public boolean getIsTimeDatum() { return timeDatum; }

    public void setDatumAttribute(String type)
    {
        int index = getStringIndex(type, dAttributes);
        if(index == dAttributeIndex)
            return;

        dAttributeIndex = index;

        reconfigureUI();
	}

    public void reconfigureUI()
    {
        // If Datum Index is NONE disable velocity bean
        if(dAttributeIndex <= NONE)
        {
            velocitySelectBean.setEditable(false);
            velocityBean.setEditable(false);
            setVelocity(0.0f);
            datumBean.setEditable(false);
            setDatum(0.0f);
        }
        else if(dAttributeIndex == USER_SPECIFIED) // User Specified
        {
            velocitySelectBean.setEditable(true);
            velocityBean.setEditable(true);
            datumBean.setEditable(true);
        }
        else
        {
            velocitySelectBean.setEditable(true);
            datumBean.setEditable(false);
        }
        if(timeDatum)
        {
            velocitySelectBean.setEditable(false);
            velocityBean.setEditable(false);
        }
    }

    public void setVelocityAttribute(String type)
    {
        int index = getStringIndex(type, vAttributes);
        if(index == vAttributeIndex)
            return;
        if(index == 0) // User Specified
            velocityBean.setEditable(true);
        else
        {
            setVelocity(0.0f);
            velocityBean.setEditable(false);
        }

        vAttributeIndex = index;
    }

    public void setDatum(double datum)
    {
        this.datum = datum;
    }

    public void setVelocity(double velocity)
    {
        this.velocity = velocity;
    }

    public String getDatumAttribute()
    {
        if(dAttributes == null)
            return attributeStrings[0];
        else
            return dAttributes[dAttributeIndex];
    }

    public String getVelocityAttribute()
    {
        if(vAttributes == null)
            return attributeStrings[0];
        else
            return vAttributes[vAttributeIndex];
	}

    public double getDatum() { return datum; }
    public double getVelocity() { return velocity; }

    static int getStringIndex(String string, String[] strings)
    {
        for(int n = 0; n < strings.length; n++)
            if(strings[n] == string)return n;
        return -1;
    }
}
