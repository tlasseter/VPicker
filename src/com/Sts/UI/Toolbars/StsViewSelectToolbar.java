package com.Sts.UI.Toolbars;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.io.*;

public class StsViewSelectToolbar extends StsToolbar implements StsSerializable
{
    transient private StsWin3dBase win3d;
    transient private StsModel model;
//    transient private StsGLPanel3d glPanel3d;
//    transient private StsWindowActionManager actionManager;

    // transient private StsViewItem viewItem = null;
    transient private StsViewItem[] singleViewItems = null;
    transient private StsViewItem[] selectedItems;
    transient private StsButton buttonChangeView, buttonAuxView, buttonFullView, buttonTimeSeriesView, buttonGaugeView;
    transient StsViewItem viewItemPlot;

    transient StsStringFieldBean viewNameBean;

    public static final String NAME = "Select View Toolbar";
    public static final boolean defaultFloatable = true;

    /** button filenames (also used as unique identifier button names) */
    public static final String RESET_VIEW = "ResetView";
    public static final String VIEW_CHANGE = "changeView";
    public static final String VIEW_AUX = "auxView";
    public static final String VIEW_FULL = "fullView";
    public static final String VIEW_TIME_SERIES = "timeSeriesView";
    public static final String VIEW_GAUGES = "gaugeView";

    public StsViewSelectToolbar()
    {
        super(NAME);
    }

    public StsViewSelectToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        this.win3d = win3d;
        this.model = win3d.getModel();
//        this.glPanel3d = win3d.getSingleViewGlPanel3d();
//        actionManager = glPanel3d.getActionManager();

        add(new StsButton(RESET_VIEW, "Reset Current View", win3d, "setDefaultView"));

        buttonChangeView = new StsButton(VIEW_CHANGE, "Change the view in this window.", this, "changeViewPopup");
        add(buttonChangeView);

        addSeparator();

        buttonAuxView = new StsButton(VIEW_AUX, "Copy view to new window.", this, "auxViewPopup");
        buttonFullView = new StsButton(VIEW_FULL, "Copy fully functional view to new window.", this, "fullViewPopup");
        add(buttonAuxView);
        add(buttonFullView);

        addSeparator();

        buttonTimeSeriesView = new StsButton(VIEW_TIME_SERIES, "Create a time series plot.", this, "timeSeriesView");
        add(buttonTimeSeriesView);
        buttonGaugeView = new StsButton(VIEW_GAUGES, "Create a gauge window.", this, "gaugeView");
        add(buttonGaugeView);

        addCloseIcon(win3d);

        setMinimumSize();

        //       setViewMethod = StsMethod.getInstanceMethod(StsWin3dBase.class, "checkAddView", new Class[] { Class.class, String.class } );
        return true;
    }

    public void changeViewPopup()
    {
        viewPopup(buttonChangeView, "Change View", "changeView");
    }

    public void auxViewPopup()
    {
        viewPopup(buttonAuxView, "New Aux View", "auxView");
    }

    public void fullViewPopup()
    {
        viewPopup(buttonFullView, "New Full View", "fullView");
    }

    private void viewPopup(StsButton button, String title, String methodName)
    {
        JPopupMenu popupMenu = buildViewPopupMenu(title, methodName);
        Point location = button.getLocation();
        popupMenu.show(button, location.x, location.y);
    }

    public void changeView(StsViewItem viewItem)
    {
        win3d.changeView(viewItem);
    }

    public void auxView(StsViewItem viewItem)
    {
        model.createAuxWindow(win3d, viewItem);
    }

    public void fullView(StsViewItem viewItem)
    {
        model.createFullWindow(viewItem);
    }

    public void timeSeriesView()
    {
       win3d.createTimeSeriesView();
    }

    public void gaugeView()
    {
       win3d.createGaugeView();
    }

    class ViewSelectDialog extends StsSimpleOkCancelDialog
    {
        StsJPanel viewItemsPanel = new StsJPanel();
        StsBooleanFieldBean[] viewSelectBeans;

        private String viewName;

        public ViewSelectDialog()
        {
            super(win3d, "SelectView");
        }

        public void layoutPanel(StsJPanel panel)
        {
            constructViewList();
            int nViewItems = singleViewItems.length;
            viewSelectBeans = new StsBooleanFieldBean[nViewItems];
            for (int n = 0; n < nViewItems; n++)
            {
                StsViewItem viewItem = singleViewItems[n];
                ViewType viewType = new ViewType(viewItem);
                viewSelectBeans[n] = viewType.booleanFieldBean;
            }
            panel.gbc.anchor = GridBagConstraints.WEST;
            panel.add(viewSelectBeans);
            viewNameBean = new StsStringFieldBean(this, "viewName", "View name");
            panel.add(viewNameBean);
        }

        public StsViewItem[] getSelectedViewItems()
        {
            return selectedItems;
        }

        public String getViewName()
        {
            return viewName;
        }

        public void setViewName(String viewName)
        {
            this.viewName = viewName;
        }
    }

    private JPopupMenu buildViewPopupMenu(String title, String methodName)
    {
        JPopupMenu popupMenu = new JPopupMenu(title);

        // create viewItems for single view types
        constructViewList();
        int nSingleViewItems = singleViewItems.length;
        StsMultiViewType[] multiViewTypes = getMultiViewTypes();
        int nMultiViewItems = multiViewTypes.length;
        StsViewItem[] multiViewItems = new StsViewItem[nMultiViewItems];

        for (int n = 0; n < nMultiViewItems; n++)
        {
            StsMultiViewType multiViewType = multiViewTypes[n];
            String name = multiViewType.name;
            StsViewItem[] viewItems = multiViewType.viewItems;
            String[] viewNames = StsViewItem.getViewNames(viewItems);
            Class[] viewClasses = StsViewItem.getViewClasses(viewItems);
            multiViewItems[n] = new StsViewItem(name, viewClasses, viewNames);
        }

        JMenuItem menuLabel = new JMenuItem("Select/create view below");
        menuLabel.setFont(new Font("Dialog", Font.BOLD, 12));
        popupMenu.add(menuLabel);
        popupMenu.addSeparator();

        StsMenuItem viewMenuItem;
        for (int n = 0; n < nSingleViewItems; n++)
        {
            StsViewItem viewItem = singleViewItems[n];
            viewMenuItem = new StsMenuItem();
            viewMenuItem.setMenuActionListener(viewItem.name, this, methodName, viewItem);
            popupMenu.add(viewMenuItem);
        }

        if(nMultiViewItems > 0)
        {
            popupMenu.addSeparator();

            for (int n = 0; n < nMultiViewItems; n++)
            {
                StsViewItem viewItem = multiViewItems[n];
                viewMenuItem = new StsMenuItem();
                viewMenuItem.setMenuActionListener(viewItem.name, this, methodName, viewItem);
                popupMenu.add(viewMenuItem);
            }
        }
        popupMenu.addSeparator();

        viewMenuItem = new StsMenuItem();
        viewMenuItem.setMenuActionListener("New multiView...", this, "newMultiView", methodName);
        popupMenu.add(viewMenuItem);

        return popupMenu;
    }

    public void newMultiView(String methodName)
    {
        ViewSelectDialog viewSelectDialog = new ViewSelectDialog();
        viewSelectDialog.setVisible(true);
        StsViewItem[] viewItems = selectedItems;
        selectedItems = null;
        if (viewSelectDialog.wasCanceled()) return;
        String viewName = viewSelectDialog.viewName;
        Class[] viewClasses = StsViewItem.getViewClasses(viewItems);
        String[] viewNames = StsViewItem.getViewNames(viewItems);
        StsViewItem multiViewItem = new StsViewItem(viewName, viewClasses, viewNames);
        StsMethod method = new StsMethod(this, methodName, multiViewItem);
        method.invokeInstanceMethod();
        StsMultiViewType[] multiViewTypes = getMultiViewTypes();
        StsMultiViewType newMultiViewType = new StsMultiViewType(viewItems, viewName);
        multiViewTypes = (StsMultiViewType[])StsMath.arrayAddElement(multiViewTypes, newMultiViewType);
        MultiViewTypeSet multiViewTypeSet = new MultiViewTypeSet(multiViewTypes);
        setMultiViewTypes(multiViewTypeSet);
    }

    public StsMultiViewType[] getMultiViewTypes()
    {
        try
        {
            String pathname = getMultiViewTypeSetsPathname();
            MultiViewTypeSet multiViewTypeSet = new MultiViewTypeSet();
            if (!StsParameterFile.initialReadObjectFields(pathname, multiViewTypeSet)) return new StsMultiViewType[0];
            StsMultiViewType[] viewTypes = multiViewTypeSet.viewTypes;
            if(viewTypes == null) viewTypes = new StsMultiViewType[0];
            for(StsMultiViewType viewType : viewTypes)
                viewType.initializeViewClasses();
            return viewTypes;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getViewTypes", e);
            return new StsMultiViewType[0];
        }
    }

    public void setMultiViewTypes(MultiViewTypeSet multiViewTypeSet)
    {
        try
        {
            String pathname = getMultiViewTypeSetsPathname();
 //           StsToolkit.writeObjectXML(multiViewTypeSet, pathname);
            StsParameterFile.writeObjectFields(pathname, multiViewTypeSet);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getViewTypes", e);
        }
    }

    private String getMultiViewTypeSetsPathname()
    {
        return System.getProperty("user.home") + File.separator + "S2SCache" + File.separator + "s2s.viewTypes";
    }

    class ViewType
    {
        StsBooleanFieldBean booleanFieldBean;
        boolean selected = false;
        StsViewItem viewItem;

        ViewType(StsViewItem viewItem)
        {
            this.viewItem = viewItem;
            booleanFieldBean = new StsBooleanFieldBean(this, "selected", viewItem.name);
        }

        public boolean getSelected()
        {
            return selected;
        }

        public void setSelected(boolean selected)
        {
            this.selected = selected;
            if(selected)
                selectedItems = (StsViewItem[])StsMath.arrayAddElement(selectedItems, viewItem);
            else
                selectedItems = (StsViewItem[])StsMath.arrayDeleteElement(selectedItems, viewItem);

            StringBuffer stringBuffer = new StringBuffer();
            int nItems = selectedItems.length;
            for(int n = 0; n < nItems-1; n++)
                stringBuffer.append(selectedItems[n].shortName + "-");
            stringBuffer.append(selectedItems[nItems-1].shortName);
            String viewName = stringBuffer.toString();
            viewNameBean.setValue(viewName);
            viewNameBean.setValueInPanelObject();
        }
    }

    class MultiViewTypeSet implements Serializable
    {
        StsMultiViewType[] viewTypes;

        public MultiViewTypeSet()
        {
        }

        public MultiViewTypeSet(StsMultiViewType[] viewTypes)
        {
            this.viewTypes = viewTypes;
        }

    }

    public boolean selectToolbarItem(String itemName, boolean selected)
    {
        return false;
    }

    private void constructViewList()
    {
        singleViewItems = StsView.constructViewList(model);
//        StsViewItem.setSelectedMethod(setViewMethod);
    }

    /*
        public void setButtonEnabled(StsView view, boolean enabled)
        {
            if(     view instanceof StsView3d)     viewItem3d.selected();
            else if(view instanceof StsViewCursor) viewItem2d.selected();
            else if(view instanceof StsViewXP)     viewItemXP.selected();
            else if(view instanceof StsViewGather3d) viewItemGather.selected();
            else if(view instanceof StsViewSemblance3d) viewItemSemblance.selected();
        }
    */

    public StsViewItem[] getSingleViewItems()
    {
        constructViewList();
        return singleViewItems;
    }

    public boolean forViewOnly() { return true; }
     
/*
    public void setViewItem(StsViewItem viewItem)
    {
        this.viewItem = viewItem;
//        viewItem.selected();
    }

    public StsViewItem getViewItem() { return viewItem; }
*/
}