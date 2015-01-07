package com.Sts.Types;

import com.Sts.DBTypes.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 9, 2010
 * Time: 2:53:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsSliceLayer extends StsSerialize
{
    StsEdgeLoop topLoop;
    StsEdgeLoop botLoop;

    public StsSliceLayer(StsEdgeLoop topLoop, StsEdgeLoop botLoop)
    {
        this.topLoop = topLoop;
        this.botLoop = botLoop;
    }
}
