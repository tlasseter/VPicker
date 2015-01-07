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
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;


public class StsDatumProperties extends StsPanelProperties
{
    double datum = 0.0f;
    double velocity = 1000.0f;
    int dAttributeIndex = 0;
    int vAttributeIndex = 0;

    transient StsGroupBox groupBox = null;

    transient static String[] dAttributes = null;
    transient static String[] vAttributes = null;

    transient StsComboBoxFieldBean datumDomainBean = null;
    transient StsComboBoxFieldBean datumSelectBean = null;
    transient StsComboBoxFieldBean velocitySelectBean = null;
    transient StsDoubleFieldBean datumBean;
    transient StsDoubleFieldBean velocityBean;

    transient static public final byte NONE = 0;
    transient static public final byte USER_SPECIFIED = 1;
    transient static public String[] attributeStrings = new String[] {"None", "User Specified"};

    transient static public final byte TIME_DOMAIN = 0;
    transient static public final byte DEPTH_DOMAIN = 1;
    transient static public String[] domainStrings = new String[] {"Time", "Depth"};
    int datumDomain = TIME_DOMAIN;

	static private final String title = "Datum Properties";

	public StsDatumProperties()
	{
	}

	public StsDatumProperties(String fieldName)
	{
		super(title, fieldName);
	}

	public StsDatumProperties(StsObject parentObject, StsDatumProperties defaultProperties, String fieldName)
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

        String[] velocityAttributes;
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
           if(datumDomain == TIME_DOMAIN)
               datumAttributes = ((StsPreStackLineSet) parentObject).lines[0].getTimeAttributes();
           else
               datumAttributes = ((StsPreStackLineSet) parentObject).lines[0].getDistanceAttributes();
        }
        else if(parentObject instanceof StsVsp)
        {
            if(datumDomain == TIME_DOMAIN)
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
        datumDomainBean = new StsComboBoxFieldBean(this, "datumDomain", "Datum Domain", domainStrings);
        datumSelectBean = new StsComboBoxFieldBean(this, "datumAttribute", "Correct Traces to Datum:", dAttributes);
        velocitySelectBean = new StsComboBoxFieldBean(this, "velocityAttribute", "Velocity for Datum Correction:", vAttributes);
        propertyBeans = new StsFieldBean[]
        {
           datumDomainBean,
           datumSelectBean,
           datumBean = new StsDoubleFieldBean(this, "datum", -10000, 10000, "Datum:", false),
           velocitySelectBean,
           velocityBean = new StsDoubleFieldBean(this, "velocity", -10000, 10000, "Velocity (m-ft/sec):", false)
        };
        reconfigureUI();
    }

    public int getDatumDomainIndex() {return datumDomain;}
    public String getDatumDomain() { return domainStrings[datumDomain];}
    public void setDatumDomain(String type)
    {
        int index = getStringIndex( type, domainStrings);
        if(index == datumDomain)
            return;
        datumDomain = index;

        getDatumAttributes();
        datumSelectBean.setListItems(dAttributes);
        if( dAttributeIndex > 1)
            dAttributeIndex = 0;
        datumSelectBean.setSelectedIndex(dAttributeIndex);
        velocitySelectBean.setSelectedIndex(vAttributeIndex);
        reconfigureUI();
    }

    public String getDatumAttribute()
    {
        if(dAttributes == null)
            return attributeStrings[0];
        else
            return dAttributes[dAttributeIndex];
    }

    public void setDatumAttribute(String type)
    {
        int index = getStringIndex(type, dAttributes);
        if(index == dAttributeIndex)
            return;

        dAttributeIndex = index;
        reconfigureUI();
	}

    public String getVelocityAttribute()
    {
        if(vAttributes == null)
            return attributeStrings[0];
        else
            return vAttributes[vAttributeIndex];
	}

    public void setVelocityAttribute(String type)
    {
        int index = getStringIndex(type, vAttributes);
        if(index == vAttributeIndex)
            return;

        vAttributeIndex = index;
        reconfigureUI();
    }

    public double getDatum() { return datum; }
    public void setDatum(double datum) {this.datum = datum;}

    public double getVelocity() { return velocity; }
    public void setVelocity(double velocity) {this.velocity = velocity;}

    public void reconfigureUI()
    {
        if( dAttributeIndex <= NONE)
        {
            velocitySelectBean.setEditable(false);
            velocityBean.setEditable(false);
            datumBean.setEditable(false);
            return;
        }

        if(dAttributeIndex == USER_SPECIFIED)
            datumBean.setEditable(true);
        else
            datumBean.setEditable(false);

        if(datumDomain == TIME_DOMAIN)
        {
            if( !isDepth)
            {
                velocitySelectBean.setEditable(false);
                velocityBean.setEditable(false);
            }
            else
            {
                velocitySelectBean.setEditable(true);
                if( vAttributeIndex == 0)
                    velocityBean.setEditable(true);
                else
                    velocityBean.setEditable(false);
            }
        }
        else
        {
            if( !isDepth)
            {
                velocitySelectBean.setEditable(true);
                if( vAttributeIndex == 0)
                    velocityBean.setEditable(true);
                else
                    velocityBean.setEditable(false);
            }
            else
            {
                velocitySelectBean.setEditable(false);
                velocityBean.setEditable(false);
            }
        }
    }

    static int getStringIndex(String string, String[] strings)
    {
        for(int n = 0; n < strings.length; n++)
            if(strings[n] == string)return n;
        return -1;
    }
}
