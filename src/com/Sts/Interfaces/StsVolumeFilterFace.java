package com.Sts.Interfaces;

import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 13, 2007
 * Time: 7:23:00 PM
 * To change this template use File | Settings | File Templates.
 */
public interface StsVolumeFilterFace
{
    public byte[] processBytePlaneData(int dir, int nCenterPlane, StsSeismicVolume volume, StsVirtualVolume virtualVolume);
    public float[] processFloatPlaneData(int dir, int nCenterPlane, StsSeismicVolume volume, StsVirtualVolume virtualVolume);
    public StsJPanel getFilterPanel();
}
