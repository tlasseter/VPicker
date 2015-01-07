package com.Sts.Types;

import com.Sts.DBTypes.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 11, 2009
 * Time: 11:55:17 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsLogTimeVector extends StsLogVector
{
    private long time;

    public StsLogTimeVector(long time, float[] values)
    {
        super(values);
        this.time = time;
    }

    public long getTime()
    {
        return time;
    }

    public void setTime(long time)
    {
        this.time = time;
    }
}
