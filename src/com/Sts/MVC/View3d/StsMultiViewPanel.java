package com.Sts.MVC.View3d;

import com.Sts.DB.*;
import com.Sts.UI.MultiSplitPane.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 19, 2008
 * Time: 9:15:55 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsMultiViewPanel extends StsMultiSplitPane implements StsSerializable
{
    String viewName;
    String[] panelNames;
    transient StsWin3dBase win3dBase;
    transient StsMultiSplitPaneLayout layout;
    static public final Dimension size00 = new Dimension(0, 0);

    public StsMultiViewPanel()
    {
    }

    private StsMultiViewPanel(String viewName, StsWin3dBase win3dBase)
    {
        this.viewName = viewName;
        this.win3dBase = win3dBase;
        constructPanelNames();
        layoutPanel(panelNames);
        addComponents(win3dBase.glPanel3ds, panelNames);
    }

    private void constructPanelNames()
    {
        int nViews = win3dBase.glPanel3ds.length;
        panelNames = new String[nViews];
        for (int n = 0; n < nViews; n++)
            panelNames[n] = "View-" + n;
    }

    protected void layoutPanel(String[] viewNames)
    {
        layout = this.getMultiSplitLayout();

        int nViews = viewNames.length;
        if (nViews == 1)
            layout.addRootLeaf(viewNames[0]);
        else
            layout.addEqualRootRowSplit(viewNames);
        setMinimumSize(size00);
    }

    static public StsMultiViewPanel constructor(String viewName, StsWin3dBase win3dBase)
    {
        try
        {
            return new StsMultiViewPanel(viewName, win3dBase);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsMultiViewPanel.class, "constructor", e);
            return null;
        }
    }

    public String getViewName() { return viewName; }

    public void addComponents(StsGLPanel3d[] glPanel3ds, String[] viewNames)
    {
        for (int n = 0; n < glPanel3ds.length; n++)
            add(glPanel3ds[n], viewNames[n]);
    }
/*
    public Dimension size()
    {
        Dimension size = super.size();
        StsException.systemDebug(this, "size", "width: " + size.width + " height: " + size.height);
        return size;
    }

    public void setSize(Dimension d)
    {
        super.resize(d);
        StsException.systemDebug(this, "setSize", "width: " + d.width + " height: " + d.height);
    }

    public void reshape(int x, int y, int width, int height)
    {
        super.reshape(x, y, width, height);
        StsException.systemDebug(this, "reshape", "width: " + width + " height: " + height);
    }
*/
}