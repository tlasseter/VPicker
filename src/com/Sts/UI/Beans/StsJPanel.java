package com.Sts.UI.Beans;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * Extends JPanel with GridBagLayout and ability to lay out row and column elements easily.
 * These are the default GridBagConstraints:
 * gridx = RELATIVE;
 * gridy = RELATIVE;
 * gridwidth = 1;
 * gridheight = 1;
 * <p/>
 * weightx = 1;
 * weighty = 1;
 * anchor = CENTER;
 * fill = BOTH;
 * <p/>
 * insets = new Insets(0, 0, 0, 0);
 * ipadx = 0;
 * ipady = 0;
 */

public class StsJPanel extends JPanel
{
    public GridBagConstraints gbc;

    {
        gbc = new GridBagConstraints();
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.fill = GridBagConstraints.BOTH;
    }

    private Component vertStrut10 = Box.createVerticalStrut(10);
    //	private Component horizGlue = Box.createHorizontalGlue();
    public boolean addVerticalSpacer = false;
    public StsFieldBean[] fieldBeans;
    protected Object[] panelObjects = null; // Object from which this panel get/set(s) values

    static public final int HORIZONTAL = GridBagConstraints.HORIZONTAL;
    static public final int VERTICAL = GridBagConstraints.VERTICAL;
    static public final int BOTH = GridBagConstraints.BOTH;
    static public final int NONE = GridBagConstraints.NONE;

    static public final int CENTER = GridBagConstraints.CENTER;
    static public final int NORTH = GridBagConstraints.NORTH;
    static public final int SOUTH = GridBagConstraints.SOUTH;
    static public final int EAST = GridBagConstraints.EAST;
    static public final int WEST = GridBagConstraints.WEST;
    static public final int LINE_END = GridBagConstraints.LINE_END;
    static public final int LINE_START = GridBagConstraints.LINE_START;

    //	public boolean isOnObjectPanel = false;  // indicates this is on objectPanel and should be laid out accordingly
    public StsJPanel()
    {
        initializeLayout(false);
    }

    public StsJPanel(boolean addInsets)
    {
        initializeLayout(addInsets);
    }
	public StsJPanel(LayoutManager l)
	{
	 setLayout(l);
    }

    public static StsJPanel noInsets() {return new StsJPanel(false);}

    public static StsJPanel addInsets() {return new StsJPanel(true);}

    void initializeLayout(boolean addInsets)
    {
        try
        {
            setLayout(new GridBagLayout());
            reinitializeLayout();
            if (addInsets)
                gbc.insets = new Insets(4, 2, 2, 4);
        }
        catch (Exception e)
        {
            StsException.systemError("StsJPanel.initializeLayout() failed.");
        }
    }

    public void reinitializeLayout()
    {
        gbc.gridx = 0;
        gbc.gridy = 0;
        gbc.gridwidth = 1;
        gbc.weightx = 1.0;
        gbc.weighty = 1.0;
        gbc.anchor = CENTER;
        gbc.fill = NONE;
        //		gbc.fill = GridBagConstraints.HORIZONTAL;
    }

    /*
    public void addFieldBean(StsFieldBean fieldBean)
    {
    if(fieldBean == null) return;

    Component[] components = fieldBean.getBeanComponents();
    components = removeNullComponents(components);
    addComponents(components);
    }
    */

    /*
    public void addFieldBeanPanel(StsFieldBean fieldBean)
    {
    if(fieldBean == null) return;

    JPanel panel = fieldBean.getPanel();
    gbc.gridwidth = 1;
    gbc.gridx = 0;
    super.add(panel, gbc);
    gbc.gridy += 1;
    if(addVerticalSpacer)
    {
    super.add(vertStrut10, gbc);
    gbc.gridy += 1;
    }
    }
    */

    /** Use this routine to add a series of components in a single vertical column */
    public Component add(Component component)
    {
        checkAddFieldBeans(component);
        gbc.gridx = 0;
        super.add(component, gbc);
        gbc.gridy += 1;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;
        return component;
    }

    public void addEmptyRows(int nRows)
    {
        gbc.gridy += nRows;
    }

    private void checkAddFieldBean(StsFieldBean fieldBean)
    {
        fieldBeans = (StsFieldBean[])StsMath.arrayAddElement(fieldBeans, fieldBean, StsFieldBean.class);
    }

    private void checkAddFieldBeans(Component component)
    {
        if (!(component instanceof StsJPanel)) return;
        fieldBeans = (StsFieldBean[]) StsMath.arrayAddArray(fieldBeans, ((StsJPanel) component).fieldBeans, StsFieldBean.class);
    }

    /** Use this routine to add a series of components in a single vertical column */
    public Component add(Component component, int ySpan, double weighty)
    {
        checkAddFieldBeans(component);
        int saveYSpan = gbc.gridwidth;
        double saveWeight = gbc.weighty;

        gbc.gridheight = ySpan;
        gbc.weighty = weighty;
        super.add(component, gbc);
        gbc.gridy += ySpan;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;

        gbc.gridheight = saveYSpan;
        gbc.weighty = saveWeight;

        return component;
    }

    public Component add(StsFieldBean fieldBean)
    {
        //		checkAddFieldBean(fieldBean);
        gbc.gridx = 0;
        return addEndRow(fieldBean);
    }

    public void add(StsFieldBean[] beans)
    {
        if (beans == null) return;
        for (int n = 0; n < beans.length; n++)
          if (beans[n] != null)
                add(beans[n]);
    }

    public void addRowOfBeans(StsFieldBean[] beans)
    {
        if (beans == null) return;
        for (int n = 0; n < beans.length - 1; n++)
            addToRow(beans[n]);
        addEndRow(beans[beans.length - 1]);
    }

    public Component addBeanPanel(StsFieldBean fieldBean)
    {
        checkAddFieldBean(fieldBean);

        int savewidth = gbc.gridwidth;
        gbc.gridwidth = fieldBean.gridwidth;
        add(fieldBean, gbc);
        gbc.gridwidth = savewidth;

        gbc.gridy += 1;
        gbc.weightx = 0.;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;
        return fieldBean;
    }

    public Component addBeanPanel(StsFieldBean fieldBean, GridBagConstraints gbc)
    {
        checkAddFieldBean(fieldBean);

        addBean(fieldBean, gbc);

        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        return fieldBean;
    }

    private void addBean(StsFieldBean fieldBean, GridBagConstraints gbc)
    {
        //gbc.gridwidth = fieldBean.gridwidth;
        add(fieldBean, gbc);
    }

    /** Adds individual components of bean to panel */
    private void addBeanComponents(StsFieldBean fieldBean)
    {
        fieldBean.addComponentsToPanel(this);
    }

    /**
     * Add the bean (an StsJPanel) as a single component to to the container.
     * If the current gridwidth is > bean.gridwidth, use it instead.
     * Increment gridx by the bean width (usually 1).
     */
    public void addBeanToRow(StsFieldBean fieldBean)
    {
        checkAddFieldBean(fieldBean);
        int width = gbc.gridwidth;  // save width so we can restore it
        gbc.gridwidth = Math.max(fieldBean.gridwidth, gbc.gridwidth);
        gbc.gridwidth = fieldBean.gridwidth;
        add(fieldBean, gbc);
        gbc.gridx += fieldBean.gridwidth;
        gbc.gridwidth = width;
    }

    /**
     * Add the bean (an StsJPanel) as a single component to to the container.
     * If the current gridwidth is > bean.gridwidth, use it instead.
     * Increment gridx by the bean width (usually 1).
     */
    public void addBeanEndRow(StsFieldBean fieldBean)
    {
        int width = gbc.gridwidth;  // save width so we can restore it
        if (fieldBean == null)
        {
            System.out.println("attempt to add null bean");
            return;
        }

        gbc.gridwidth = Math.max(fieldBean.gridwidth, gbc.gridwidth);
        checkAddFieldBean(fieldBean);
        add(fieldBean, gbc);
        gbc.gridwidth = width; // restore width
        gbc.gridx = 0;
        gbc.gridy += 1;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
    }

    /** Use this routine to add a component to the end of the current row  and start a new row.  Use addToRow if there is no next row. */
    public Component addEndRow(Component component)
    {
        checkAddFieldBeans(component);
        super.add(component, gbc);
        gbc.gridy += 1;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;
        return component;
    }

    /** Use this routine to add a component to the current row */
    public Component addEndRow(Component component, int xSpan, double weight)
    {
        checkAddFieldBeans(component);
        int saveXSpan = gbc.gridwidth;
        double saveWeight = gbc.weightx;

        gbc.gridwidth = xSpan;
        gbc.weightx = weight;
        super.add(component, gbc);
        gbc.gridy += 1;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;

        gbc.gridwidth = saveXSpan;
        gbc.weightx = saveWeight;

        return component;
    }

    /** Use this routine to add a component to the end of the current row  and start a new row.  Use addToRow if there is no next row. */
    public Component addEndRow(StsFieldBean fieldBean)
    {
        addBeanEndRow(fieldBean);
        return fieldBean;
    }

    /*
    public Component addEndRow(StsFieldBean fieldBean)
    {
        checkAddFieldBean(fieldBean);

        int savewidth = gbc.gridwidth;
        gbc.gridwidth = fieldBean.gridwidth;
        add(fieldBean, gbc);
        gbc.gridwidth = savewidth;

        gbc.gridy += 1;
        if (addVerticalSpacer)
        {
            super.add(vertStrut10, gbc);
            gbc.gridy += 1;
        }
        gbc.gridx = 0;
        return fieldBean;
    }
    */
    private void checkAddFieldBean(Component component)
    {
        if (!(component instanceof StsFieldBean)) return;
        StsFieldBean bean = (StsFieldBean) component;
        fieldBeans = (StsFieldBean[]) StsMath.arrayAddElement(fieldBeans, bean, StsFieldBean.class);
    }

    public StsFieldBean[] getFieldBeans()
    {
        return fieldBeans;
    }

    /** Use this routine to add a component to the current row */
    public Component addToRow(Component component)
    {
        checkAddFieldBeans(component);
        //		gbc.gridwidth = 1;
        super.add(component, gbc);
        gbc.gridx += gbc.gridwidth;
        return component;
    }

    /** Use this routine to add a component to the current row */
    public Component addToRow(Component component, int xSpan, double weight)
    {
        checkAddFieldBeans(component);
        int saveXSpan = gbc.gridwidth;
        double saveWeight = gbc.weightx;

        gbc.gridwidth = xSpan;
        gbc.weightx = weight;
        super.add(component, gbc);
        gbc.gridx += xSpan;

        gbc.gridwidth = saveXSpan;
        gbc.weightx = saveWeight;

        return component;
    }

    public Component addToRow(StsFieldBean fieldBean)
    {
        addBeanToRow(fieldBean);
        return fieldBean;
    }
    /*
    public Component addToRow(StsFieldBean fieldBean)
    {
        checkAddFieldBean(fieldBean);
        // addBean(fieldBean);
        addBeanComponents(fieldBean);
        return fieldBean;
    }
    */
    /*
    public void addComponents(Component[] components)
    {
    int nComponents = components.length;
    if(nComponents > 0)
    {
    gbc.gridwidth = 1;
    gbc.gridx = 0;
    super.add(components[0], gbc);
    }
    if(nComponents > 1)
    {
    gbc.gridwidth = GridBagConstraints.REMAINDER;
    gbc.gridx = 1;
    super.add(components[1], gbc);
    }
    gbc.gridy += 1;
    if(addVerticalSpacer)
    {
    super.add(vertStrut10, gbc);
    gbc.gridy += 1;
    }
    }

    private Component[] removeNullComponents(Component[] components)
    {
    if(components == null) return new Component[0];
    int nNotNull = 0;
    for(int n = 0; n < components.length; n++)
    if(components[n] != null) nNotNull++;

    if(nNotNull == components.length) return components;

    Component[] newComponents = new Component[nNotNull];
    int nn = 0;
    for(int n = 0; n < components.length; n++)
    if(components[n] != null) newComponents[nn++] = components[n];
    return newComponents;
    }
    */


    public void setEditable(boolean editable)
    {
        if (fieldBeans == null) return;
        for (int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].setEditable(editable);
    }

    public void setPanelObject(Object object)
    {
        if (fieldBeans == null) return;
        panelObjects = new Object[]{object};
        for (int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].setBeanObject(object);
    }

    public void initializeBeans()
    {
        if (fieldBeans == null) return;
        for (int n = 0; n < fieldBeans.length; n++)
            fieldBeans[n].initializeBean();
    }

    public void updateBeans()
    {
        initializeBeans();
    }

    public Object getPanelObject()
    {
        if (panelObjects == null)
            return null;
        return panelObjects[0];
    }

    public void setPanelObjects(Object[] objects)
    {
        panelObjects = objects;
    }

    public Object[] getPanelObjects()
    {
        return panelObjects;
    }

    public void objectPanelLayout(Component[] components)
    {
        gbc.fill = HORIZONTAL;
        gbc.anchor = LINE_START;
        for (int n = 0; n < components.length; n++)
            addToRow(components[n]);
        gbc.anchor = CENTER;
        gbc.fill = NONE;
    }

    public void setPreferredSize(int width, int height)
    {
        setPreferredSize(new Dimension(width, height));
    }

    static public void main(String[] args)
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                final StsJPanel panel = StsJPanel.addInsets();
                panel.setPreferredSize(100, 100);
                JFrame frame = new JFrame();
                frame.add(panel, BorderLayout.CENTER);

                MouseListener mouseListener = new MouseAdapter()
                {
                    public void mousePressed(MouseEvent e)
                    {
                        System.out.println(" mouse pressed.");
                    }

                    public void mouseReleased(MouseEvent e)
                    {
                        System.out.println("mouse released.");
                    }

                    public void mouseDragged(MouseEvent e)
                    {
                        System.out.println("mouse dragged.");
                    }
                };
                ComponentListener componentListener = new ComponentAdapter()
                {
                    int width
                        ,
                        height;

                    public void componentResized(ComponentEvent e)
                    {
                        System.out.println("component resized.");
                        JPanel resizedPanel = (JPanel) e.getSource();
                        int newWidth = resizedPanel.getWidth();
                        int newHeight = resizedPanel.getHeight();
                        System.out.println("resizedWidth: " + newWidth + " resizedHeight: " + newHeight);
                        System.out.println("panelWidth: " + panel.getWidth() + " paneldHeight: " + panel.getHeight());
                        System.out.println("oldWidth: " + width + " oldHeight: " + height);
                        boolean changed = newWidth != width || newHeight != height;
                        if (changed)
                        {
                            width = newWidth;
                            height = newHeight;
                            System.out.println("size changed");
                        }
                        else
                            System.out.println("size not changed");
                    }

                    /**
                     * Invoked when the component's position changes.
                     */
                    public void componentMoved(ComponentEvent e)
                    {

                        System.out.println("component moved.");
                    }
                };
                panel.addMouseListener(mouseListener);
                panel.addComponentListener(componentListener);
                frame.setVisible(true);
            }

        };
        StsToolkit.runLaterOnEventThread(runnable);
    }
}
