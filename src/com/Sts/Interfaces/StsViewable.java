package com.Sts.Interfaces;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 28, 2008
 * Time: 8:06:21 AM
 * To change this template use File | Settings | File Templates.
 */

/** Classes which implement StsViewable require specialized views.  getViewClasses() returns an array of viewClasses
 * which can display this class of objects.  This interface is used in assembling possible views in the view dropdown list.
 */
public interface StsViewable
{
    public Class[] getViewClasses();
}
