package com.Sts.Types;

import com.Sts.DBTypes.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:14:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsNullAttribute extends StsDisplaySeismicAttribute
{
    public String getName() { return StsPreStackLineSet.ATTRIBUTE_NONE_STRING; }

    public int getIndex() { return StsPreStackLineSet.ATTRIBUTE_NONE; }

    public boolean computeBytes()
    {
        return false;
    }

    public void addColorscale() { }
}
