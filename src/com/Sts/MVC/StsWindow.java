package com.Sts.MVC;

import com.Sts.DB.*;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Oct 6, 2008
 * Time: 9:20:47 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsWindow extends JFrame implements StsSerializable
{
    transient public StsActionManager actionManager;
    public transient StsModel model;

    public StsWindow()
    {
    }

    public StsWindow(StsModel model, StsActionManager actionManager)
    {
        this.model = model;
        this.actionManager = actionManager;
    }

    public StsModel getModel()
    {
        return model;
    }

    public GraphicsDevice getGraphicsDevice()
    {
        GraphicsConfiguration gConfig = getGraphicsConfiguration();
        if(gConfig == null) return null;
        return gConfig.getDevice();
    }
}
