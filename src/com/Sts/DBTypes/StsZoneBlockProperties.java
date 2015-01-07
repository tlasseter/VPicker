package com.Sts.DBTypes;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 26, 2010
 * Time: 10:37:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsZoneBlockProperties extends StsObject
{
    StsZoneBlock zoneBlock;
    String propertyName;
    byte type;
    transient float[][][] subZoneValues;
    transient float[] constantSubZoneValues;
    transient float constantValue;

    static public final byte PROP_NONE = 0;
    static public final byte PROP_CONSTANT = 1;
    static public final byte PROP_SUBZONE_CONSTANT = 2;
    static public final byte PROP_VARIABLE = 3;

    public StsZoneBlockProperties(StsZoneBlock zoneBlock, String propertyName)
    {
        this.zoneBlock = zoneBlock;
        this.propertyName = propertyName;
    }

    public void setConstant(float value)
    {
        type = PROP_CONSTANT;
        constantValue = value;
    }

    public void setSubZoneConstant(float[] values)
    {
        type = PROP_SUBZONE_CONSTANT;
        constantSubZoneValues = values;
    }

    public void setVariable(float[][][] values)
    {
        type = PROP_VARIABLE;
        subZoneValues = values;
    }
}
