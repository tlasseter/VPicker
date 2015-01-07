package com.Sts.DBTypes;

import com.Sts.Interfaces.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 8, 2007
 * Time: 12:57:25 PM
 * To change this template use File | Settings | File Templates.
 */

/** This class contains the instances of prestack volumes which are velocity, stack, and semblance volumes.
 *  They could be included in the StsSeismicVolumeClass since they are poststack volumes,
 *  but for the moment they have been separated out.
 */
public class StsPreStackVolumeClass extends StsSeismicVolumeClass implements StsClassCursor3dTextureDisplayable
{
	public StsPreStackVolumeClass()
	{
		super();
        userName = "Pre-Stack Volume";
	}
}
