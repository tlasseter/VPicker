package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:03:44 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsDisplayPreStackAttribute extends StsDisplaySeismicAttribute
{
    protected StsPreStackLineSet lineSet;

    public StsDisplayPreStackAttribute(StsPreStackLineSet lineSet, StsModel model)
    {
        this.model = model;
        this.lineSet = lineSet;
    }
    public StsColorscale getColorscaleWithName(String name)
    {
        return lineSet.getColorscaleWithName(name);
    }

    public void addColorscale(StsColorscale colorscale)
    {
        lineSet.addColorscale(colorscale);
    }
}
