package com.Sts.Types;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 11, 2010
 * Time: 5:44:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsPropertyTypeDefinition
{
    public String eclipseName;
    public String name;
    public String spectrumName = null;

    public StsPropertyTypeDefinition(String name, String eclipseName, String spectrumName)
    {
        this.name = name;
        this.eclipseName = eclipseName;
        this.spectrumName = spectrumName;
    }

    public StsPropertyTypeDefinition(String name, String eclipseName)
    {
        this.name = name;
        this.eclipseName = eclipseName;
    }

    public String toString() { return name; }
}
