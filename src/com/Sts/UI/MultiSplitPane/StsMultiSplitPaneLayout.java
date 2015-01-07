package com.Sts.UI.MultiSplitPane;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 16, 2008
 * Time: 8:28:12 AM
 * To change this template use File | Settings | File Templates.
 */


import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.beans.*;
import java.util.*;


/**
 * The MultiSplitLayout layout manager recursively arranges its
 * components in row and column groups called "Splits".  Elements of
 * the layout are separated by gaps called "Dividers".  The overall
 * layout is defined with a simple tree model whose nodes are
 * instances of MultiSplitLayout.Split, MultiSplitLayout.Divider,
 * and MultiSplitLayout.Leaf. Named Leaf nodes represent the space
 * allocated to a component that was added with a constraint that
 * matches the Leaf's name.  Extra space is distributed
 * among row/column siblings according to their 0.0 to 1.0 weight.
 * If no weights are specified then the last sibling always gets
 * all of the extra space, or space reduction.
 * <p/>
 * <p/>
 * Although MultiSplitLayout can be used with any Container, it's
 * the default layout manager for MultiSplitPane.  MultiSplitPane
 * supports interactively dragging the Dividers, accessibility,
 * and other features associated with split panes.
 * <p/>
 * <p/>
 * All properties in this class are bound: when a properties value
 * is changed, all PropertyChangeListeners are fired.
 *
 * Based on swingx MultiSplitPane, but heavily modified
 */

public class StsMultiSplitPaneLayout implements LayoutManager
{
    private final PropertyChangeSupport pcs = new PropertyChangeSupport(this);
    private Leaf root;
    private int dividerSize;
    private boolean floatingDividers = true;
    private boolean draggingDivider = false;
    private ArrayList<Leaf> leafList = new ArrayList();

    static public byte NONE = 0;
    static public byte ROOT = 1;
    static public byte LEAF = 2;
    static public byte ROW = 3;
    static public byte COL = 4;
    static public byte DIVIDER = 5;

    static final boolean debug = false;

    public StsMultiSplitPaneLayout()
    {
//        this.root = new Leaf();
        this.dividerSize = UIManager.getInt("SplitPane.dividerSize");
        if (this.dividerSize == 0)
        {
            this.dividerSize = 7;
        }
    }

    public StsMultiSplitPaneLayout(byte splitType)
    {
        this.root = new Leaf("ROOT", splitType);
        this.dividerSize = UIManager.getInt("SplitPane.dividerSize");
        if (this.dividerSize == 0)
        {
            this.dividerSize = 7;
        }
    }

    public void addPropertyChangeListener(PropertyChangeListener listener)
    {
        if (listener != null)
        {
            pcs.addPropertyChangeListener(listener);
        }
    }

    public void removePropertyChangeListener(PropertyChangeListener listener)
    {
        if (listener != null)
        {
            pcs.removePropertyChangeListener(listener);
        }
    }

    public PropertyChangeListener[] getPropertyChangeListeners()
    {
        return pcs.getPropertyChangeListeners();
    }

    private void firePCS(String propertyName, Object oldValue, Object newValue)
    {
        if(debug) StsException.systemDebug(this, "firePCS", "floatingDividers changed from: " + oldValue + " to " + newValue);
//        if (!(oldValue != null && newValue != null && oldValue.equals(newValue)))
        {
            pcs.firePropertyChange(propertyName, oldValue, newValue);
        }
    }

    /**
     * Return the root of the tree of Split, Leaf, and Divider nodes
     * that define this layout.
     *
     * @return the value of the model property
     */
    public Node getRoot()
    { return root; }

    /**
     * Returns the width of Dividers in Split rows, and the height of
     * Dividers in Split columns.
     *
     * @return the value of the dividerSize property
     * @see #setDividerSize
     */
    public int getDividerSize()
    { return dividerSize; }

    /**
     * Sets the width of Dividers in Split rows, and the height of
     * Dividers in Split columns.  The default value of this property
     * is the same as for JSplitPane Dividers.
     *
     * @param dividerSize the size of dividers (pixels)
     * @throws IllegalArgumentException if dividerSize < 0
     * @see #getDividerSize
     */
    public void setDividerSize(int dividerSize)
    {
        if (dividerSize < 0)
        {
            throw new IllegalArgumentException("invalid dividerSize");
        }
        int oldDividerSize = this.dividerSize;
        this.dividerSize = dividerSize;
        firePCS("dividerSize", oldDividerSize, dividerSize);
    }

    public void setDraggingDivider(boolean isDragging)
    {
        if(draggingDivider == isDragging) return;
        draggingDivider = isDragging;
        firePCS("draggingDividers", !draggingDivider, draggingDivider);
    }


    /**
     * Add a component to this MultiSplitLayout.  The
     * <code>name</code> should match the name property of the Leaf
     * node that represents the bounds of <code>child</code>.  After
     * layoutContainer() recomputes the bounds of all of the nodes in
     * the model, it will set this child's bounds to the bounds of the
     * Leaf node with <code>name</code>.  Note: if a component was already
     * added with the same name, this method does not remove it from
     * its parent.
     *
     * @param name  identifies the Leaf node that defines the child's bounds
     * @param component the component to be added
     * @see #removeLayoutComponent
     */
    public void addLayoutComponent(String name, Component component)
    {
        if (name == null)
        {
            throw new IllegalArgumentException("name not specified");
        }
        for(Leaf leaf: leafList)
        {
            if(leaf.getName().equals(name))
            {
                if(leaf.isSplit())
                    System.err.println("can't add component to  node " + leaf.getName() + " which is not a leaf.");
                else
                    leaf.setComponent(component);
                return;
            }
        };
    }

    /**
     * Removes the specified component from the layout.
     *
     * @param component the component to be removed
     * @see #addLayoutComponent
     */
    public void removeLayoutComponent(Component component)
    {
        removeLayoutLeaf(component.getName());
    }

    public void removeLayoutLeaf(String name)
    {
        for(Leaf leaf : leafList)
        {
            if(leaf.getName().equals(name))
            {
                removeLeaf(leaf);
                return;
            }
        };
    }

    public void removeLeaf(Leaf leaf)
    {
        leaf.delete();
        leafList.remove(leaf);
    }

    private Dimension preferredComponentSize(Node node)
    {
        Component component = node.getComponent();
        return (component != null) ? component.getPreferredSize() : new Dimension(0, 0);

    }

    private Dimension minimumComponentSize(Node node)
    {
        Component component = node.getComponent();
        return (component != null) ? component.getMinimumSize() : new Dimension(0, 0);

    }

    private Dimension preferredNodeSize(Node node)
    {
        if (node.type == LEAF)
        {
            return preferredComponentSize(node);
        }
        else if (node.type == DIVIDER)
        {
            int dividerSize = getDividerSize();
            return new Dimension(dividerSize, dividerSize);
        }
        else // ROW or COL split
        {
            ArrayList<Node> splitChildren = ((Leaf)node).getChildren();
            int width = 0;
            int height = 0;
            if (node.type == ROW)
            {
                for (Node splitChild : splitChildren)
                {
                    Dimension size = preferredNodeSize(splitChild);
                    width += size.width;
                    height = Math.max(height, size.height);
                }
            }
            else
            {
                for (Node splitChild : splitChildren)
                {
                    Dimension size = preferredNodeSize(splitChild);
                    width = Math.max(width, size.width);
                    height += size.height;
                }
            }
            return new Dimension(width, height);
        }
    }

    private Dimension minimumNodeSize(Node node)
    {
        if (node.type == LEAF)
        {
            Component child = node.getComponent();
            return (child != null) ? child.getMinimumSize() : new Dimension(0, 0);
        }
        else if (node.type == DIVIDER)
        {
            int dividerSize = getDividerSize();
            return new Dimension(dividerSize, dividerSize);
        }
        else
        {
            Leaf split = (Leaf)node;
            java.util.List<Node> splitChildren = split.getChildren();
            int width = 0;
            int height = 0;
            if (split.type == ROW)
            {
                for (Node splitChild : splitChildren)
                {
                    Dimension size = minimumNodeSize(splitChild);
                    width += size.width;
                    height = Math.max(height, size.height);
                }
            }
            else
            {
                for (Node splitChild : splitChildren)
                {
                    Dimension size = minimumNodeSize(splitChild);
                    width = Math.max(width, size.width);
                    height += size.height;
                }
            }
            return new Dimension(width, height);
        }
    }

    private Dimension sizeWithInsets(Container parent, Dimension size)
    {
        Insets insets = parent.getInsets();
        int width = size.width + insets.left + insets.right;
        int height = size.height + insets.top + insets.bottom;
        return new Dimension(width, height);
    }

    public Dimension preferredLayoutSize(Container parent)
    {
        Dimension size = preferredNodeSize(getRoot());
        return sizeWithInsets(parent, size);
    }

    public Dimension minimumLayoutSize(Container parent)
    {
        Dimension size = minimumNodeSize(getRoot());
        return sizeWithInsets(parent, size);
    }


    private Rectangle boundsWithYandHeight(Rectangle bounds, double y, double height)
    {
        Rectangle r = new Rectangle();
        r.setBounds((int) (bounds.getX()), (int) y, (int) (bounds.getWidth()), (int) height);
        return r;
    }

    private Rectangle boundsWithXandWidth(Rectangle bounds, double x, double width)
    {
        Rectangle r = new Rectangle();
        r.setBounds((int) x, (int) (bounds.getY()), (int) width, (int) (bounds.getHeight()));
        return r;
    }


    private void minimizeSplitBounds(Leaf split, Rectangle bounds)
    {
        Rectangle splitBounds = new Rectangle(bounds.x, bounds.y, 0, 0);
        java.util.List<Node> splitChildren = split.getChildren();
        Node lastChild = splitChildren.get(splitChildren.size() - 1);
        Rectangle lastChildBounds = lastChild.getBounds();
        if (split.isRowLayout())
        {
            int lastChildMaxX = lastChildBounds.x + lastChildBounds.width;
            splitBounds.add(lastChildMaxX, bounds.y + bounds.height);
        }
        else
        {
            int lastChildMaxY = lastChildBounds.y + lastChildBounds.height;
            splitBounds.add(bounds.x + bounds.width, lastChildMaxY);
        }
        split.setBounds(splitBounds);
    }


    private void layoutShrink(Leaf split, Rectangle bounds)
    {
        Rectangle splitBounds = split.getBounds();
        ListIterator<Node> splitChildren = split.getChildren().listIterator();
        Node lastWeightedChild = split.lastWeightedChild();

        if (split.isRowLayout())
        {
            int totalWidth = 0;          // sum of the children's widths
            int minWeightedWidth = 0;    // sum of the weighted childrens' min widths
            int totalWeightedWidth = 0;  // sum of the weighted childrens' widths
            for (Node splitChild : split.getChildren())
            {
                int nodeWidth = splitChild.getBounds().width;
                int nodeMinWidth = Math.min(nodeWidth, minimumNodeSize(splitChild).width);
                totalWidth += nodeWidth;
                if (splitChild.getWeight() > 0.0)
                {
                    minWeightedWidth += nodeMinWidth;
                    totalWeightedWidth += nodeWidth;
                }
            }

            double x = bounds.getX();
            double extraWidth = splitBounds.getWidth() - bounds.getWidth();
            double availableWidth = extraWidth;
            boolean onlyShrinkWeightedComponents =
                (totalWeightedWidth - minWeightedWidth) > extraWidth;

            while (splitChildren.hasNext())
            {
                Node splitChild = splitChildren.next();
                Rectangle splitChildBounds = splitChild.getBounds();
                double minSplitChildWidth = minimumNodeSize(splitChild).getWidth();
                double splitChildWeight = (onlyShrinkWeightedComponents)
                    ? splitChild.getWeight()
                    : (splitChildBounds.getWidth() / (double) totalWidth);

                if (!splitChildren.hasNext())
                {
                    double newWidth = Math.max(minSplitChildWidth, bounds.getMaxX() - x);
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, newWidth);
                    layout2(splitChild, newSplitChildBounds);
                }
                else if ((availableWidth > 0.0) && (splitChildWeight > 0.0))
                {
                    double allocatedWidth = Math.rint(splitChildWeight * extraWidth);
                    double oldWidth = splitChildBounds.getWidth();
                    double newWidth = Math.max(minSplitChildWidth, oldWidth - allocatedWidth);
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, newWidth);
                    layout2(splitChild, newSplitChildBounds);
                    availableWidth -= (oldWidth - splitChild.getBounds().getWidth());
                }
                else
                {
                    double existingWidth = splitChildBounds.getWidth();
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, existingWidth);
                    layout2(splitChild, newSplitChildBounds);
                }
                x = splitChild.getBounds().getMaxX();
            }
        }

        else
        {
            int totalHeight = 0;          // sum of the children's heights
            int minWeightedHeight = 0;    // sum of the weighted childrens' min heights
            int totalWeightedHeight = 0;  // sum of the weighted childrens' heights
            for (Node splitChild : split.getChildren())
            {
                int nodeHeight = splitChild.getBounds().height;
                int nodeMinHeight = Math.min(nodeHeight, minimumNodeSize(splitChild).height);
                totalHeight += nodeHeight;
                if (splitChild.getWeight() > 0.0)
                {
                    minWeightedHeight += nodeMinHeight;
                    totalWeightedHeight += nodeHeight;
                }
            }

            double y = bounds.getY();
            double extraHeight = splitBounds.getHeight() - bounds.getHeight();
            double availableHeight = extraHeight;
            boolean onlyShrinkWeightedComponents =
                (totalWeightedHeight - minWeightedHeight) > extraHeight;

            while (splitChildren.hasNext())
            {
                Node splitChild = splitChildren.next();
                Rectangle splitChildBounds = splitChild.getBounds();
                double minSplitChildHeight = minimumNodeSize(splitChild).getHeight();
                double splitChildWeight = (onlyShrinkWeightedComponents)
                    ? splitChild.getWeight()
                    : (splitChildBounds.getHeight() / (double) totalHeight);

                if (!splitChildren.hasNext())
                {
                    double oldHeight = splitChildBounds.getHeight();
                    double newHeight = Math.max(minSplitChildHeight, bounds.getMaxY() - y);
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, newHeight);
                    layout2(splitChild, newSplitChildBounds);
                    availableHeight -= (oldHeight - splitChild.getBounds().getHeight());
                }
                else if ((availableHeight > 0.0) && (splitChildWeight > 0.0))
                {
                    double allocatedHeight = Math.rint(splitChildWeight * extraHeight);
                    double oldHeight = splitChildBounds.getHeight();
                    double newHeight = Math.max(minSplitChildHeight, oldHeight - allocatedHeight);
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, newHeight);
                    layout2(splitChild, newSplitChildBounds);
                    availableHeight -= (oldHeight - splitChild.getBounds().getHeight());
                }
                else
                {
                    double existingHeight = splitChildBounds.getHeight();
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, existingHeight);
                    layout2(splitChild, newSplitChildBounds);
                }
                y = splitChild.getBounds().getMaxY();
            }
        }

        /* The bounds of the Split node root are set to be
       * big enough to contain all of its children. Since
       * Leaf children can't be reduced below their
       * (corresponding java.awt.Component) minimum sizes,
       * the size of the Split's bounds maybe be larger than
       * the bounds we were asked to fit within.
       */
        minimizeSplitBounds(split, bounds);
    }


    private void layoutGrow(Leaf split, Rectangle bounds)
    {
        Rectangle splitBounds = split.getBounds();
        ListIterator<Node> splitChildren = split.getChildren().listIterator();
        Node lastWeightedChild = split.lastWeightedChild();

        /* Layout the Split's child Nodes' along the X axis.  The bounds
       * of each child will have the same y coordinate and height as the
       * layoutGrow() bounds argument.  Extra width is allocated to the
       * to each child with a non-zero weight:
       *     newWidth = currentWidth + (extraWidth * splitChild.getWeight())
       * Any extraWidth "left over" (that's availableWidth in the loop
       * below) is given to the last child.  Note that Dividers always
       * have a weight of zero, and they're never the last child.
       */
        if (split.type == ROW)
        {
            double x = bounds.getX();
            double extraWidth = bounds.getWidth() - splitBounds.getWidth();
            double availableWidth = extraWidth;

            while (splitChildren.hasNext())
            {
                Node splitChild = splitChildren.next();
                Rectangle splitChildBounds = splitChild.getBounds();
                double splitChildWeight = splitChild.getWeight();

                if (!splitChildren.hasNext())
                {
                    double newWidth = bounds.getMaxX() - x;
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, newWidth);
                    layout2(splitChild, newSplitChildBounds);
                }
                else if ((availableWidth > 0.0) && (splitChildWeight > 0.0))
                {
                    double allocatedWidth = (splitChild.equals(lastWeightedChild))
                        ? availableWidth
                        : Math.rint(splitChildWeight * extraWidth);
                    double newWidth = splitChildBounds.getWidth() + allocatedWidth;
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, newWidth);
                    layout2(splitChild, newSplitChildBounds);
                    availableWidth -= allocatedWidth;
                }
                else
                {
                    double existingWidth = splitChildBounds.getWidth();
                    Rectangle newSplitChildBounds = boundsWithXandWidth(bounds, x, existingWidth);
                    layout2(splitChild, newSplitChildBounds);
                }
                x = splitChild.getBounds().getMaxX();
            }
        }

        /* Layout the Split's child Nodes' along the Y axis.  The bounds
       * of each child will have the same x coordinate and width as the
       * layoutGrow() bounds argument.  Extra height is allocated to the
       * to each child with a non-zero weight:
       *     newHeight = currentHeight + (extraHeight * splitChild.getWeight())
       * Any extraHeight "left over" (that's availableHeight in the loop
       * below) is given to the last child.  Note that Dividers always
       * have a weight of zero, and they're never the last child.
       */
        else
        {
            double y = bounds.getY();
            double extraHeight = bounds.getMaxY() - splitBounds.getHeight();
            double availableHeight = extraHeight;

            while (splitChildren.hasNext())
            {
                Node splitChild = splitChildren.next();
                Rectangle splitChildBounds = splitChild.getBounds();
                double splitChildWeight = splitChild.getWeight();

                if (!splitChildren.hasNext())
                {
                    double newHeight = bounds.getMaxY() - y;
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, newHeight);
                    layout2(splitChild, newSplitChildBounds);
                }
                else if ((availableHeight > 0.0) && (splitChildWeight > 0.0))
                {
                    double allocatedHeight = (splitChild.equals(lastWeightedChild))
                        ? availableHeight
                        : Math.rint(splitChildWeight * extraHeight);
                    double newHeight = splitChildBounds.getHeight() + allocatedHeight;
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, newHeight);
                    layout2(splitChild, newSplitChildBounds);
                    availableHeight -= allocatedHeight;
                }
                else
                {
                    double existingHeight = splitChildBounds.getHeight();
                    Rectangle newSplitChildBounds = boundsWithYandHeight(bounds, y, existingHeight);
                    layout2(splitChild, newSplitChildBounds);
                }
                y = splitChild.getBounds().getMaxY();
            }
        }
    }


    /* Second pass of the layout algorithm: branch to layoutGrow/Shrink
     * as needed.
     */
    private void layout2(Node node, Rectangle bounds)
    {
        if (node.type == LEAF)
        {
            Component child = node.getComponent();
            if (child != null)
            {
                child.setBounds(bounds);
            }
            node.setBounds(bounds);
        }
        else if (node.type == DIVIDER)
        {
            node.setBounds(bounds);
        }
        else // ROW or COL split
        {
            boolean grow = node.type == ROW ? (node.getBounds().width <= bounds.width)
                                              : (node.getBounds().height <= bounds.height);
            if (grow)
            {
                layoutGrow((Leaf)node, bounds);
                node.setBounds(bounds);
            }
            else
            {
                layoutShrink((Leaf)node, bounds);
                // split.setBounds() called in layoutShrink()
            }
        }
    }


    /* First pass of the layout algorithm.
     *
     * If the Dividers are "floating" then set the bounds of each
     * node to accomodate the preferred size of all of the
     * Leaf's java.awt.Components.  Otherwise, just set the bounds
     * of each Leaf/Split node so that it's to the left of (for
     * Split.isRowLayout() Split children) or directly above
     * the Divider that follows.
     *
     * This pass sets the bounds of each Node in the layout model.  It
     * does not resize any of the parent Container's
     * (java.awt.Component) children.  That's done in the second pass,
     * see layoutGrow() and layoutShrink().
     */
    private void layout1(Node node, Rectangle bounds)
    {
        if (node.type == LEAF)
        {
            node.setBounds(bounds);
        }
        else if (node.isSplit())
        {
            Leaf split = (Leaf) node;
            Iterator<Node> splitChildren = split.getChildren().iterator();
            Rectangle childBounds = null;
            int dividerSize = getDividerSize();

            /* Layout the Split's child Nodes' along the X axis.  The bounds
            * of each child will have the same y coordinate and height as the
            * layout1() bounds argument.
            *
            * Note: the column layout code - that's the "else" clause below
            * this if, is identical to the X axis (rowLayout) code below.
            */
            if (split.type == ROW)
            {
                double x = bounds.getX();
                while (splitChildren.hasNext())
                {
                    Node splitChild = splitChildren.next();
                    Divider dividerChild =
                        (splitChildren.hasNext()) ? (Divider) (splitChildren.next()) : null;

                    double childWidth = 0.0;
                    if(!draggingDivider)
                    // if (getFloatingDividers())
                    {
                        childWidth = preferredNodeSize(splitChild).getWidth();
                    }
                    else
                    {
                        if (dividerChild != null)
                        {
                            childWidth = dividerChild.getBounds().getX() - x;
                        }
                        else
                        {
                            childWidth = split.getBounds().getMaxX() - x;
                        }
                    }
                    childBounds = boundsWithXandWidth(bounds, x, childWidth);
                    layout1(splitChild, childBounds);
                    if (!draggingDivider && (dividerChild != null))
                    // if (getFloatingDividers() && (dividerChild != null))
                    {
                        double dividerX = childBounds.getMaxX();
                        Rectangle dividerBounds = boundsWithXandWidth(bounds, dividerX, dividerSize);
                        dividerChild.setBounds(dividerBounds);
                    }
                    if (dividerChild != null)
                    {
                        x = dividerChild.getBounds().getMaxX();
                    }
                }
            }

            /* Layout the Split's child Nodes' along the Y axis.  The bounds
            * of each child will have the same x coordinate and width as the
            * layout1() bounds argument.  The algorithm is identical to what's
            * explained above, for the X axis case.
            */
            else
            {
                double y = bounds.getY();
                while (splitChildren.hasNext())
                {
                    Node splitChild = splitChildren.next();
                    Divider dividerChild =
                        (splitChildren.hasNext()) ? (Divider) (splitChildren.next()) : null;

                    double childHeight = 0.0;
                    if (!draggingDivider)
                    // if (getFloatingDividers())
                    {
                        childHeight = preferredNodeSize(splitChild).getHeight();
                    }
                    else
                    {
                        if (dividerChild != null)
                        {
                            childHeight = dividerChild.getBounds().getY() - y;
                        }
                        else
                        {
                            childHeight = split.getBounds().getMaxY() - y;
                        }
                    }
                    childBounds = boundsWithYandHeight(bounds, y, childHeight);
                    layout1(splitChild, childBounds);

                    if (!draggingDivider && (dividerChild != null))
                    // if (getFloatingDividers() && (dividerChild != null))
                    {
                        double dividerY = childBounds.getMaxY();
                        Rectangle dividerBounds = boundsWithYandHeight(bounds, dividerY, dividerSize);
                        dividerChild.setBounds(dividerBounds);
                    }
                    if (dividerChild != null)
                    {
                        y = dividerChild.getBounds().getMaxY();
                    }
                }
            }
            /* The bounds of the Split node node are set to be just
            * big enough to contain all of its children, but only
            * along the axis it's allocating space on.  That's
            * X for rows, Y for columns.  The second pass of the
            * layout algorithm - see layoutShrink()/layoutGrow()
            * allocates extra space.
            */
            minimizeSplitBounds(split, bounds);
        }
    }

    /**
     * The specified Node is either the wrong type or was configured
     * incorrectly.
     */
    public static class InvalidLayoutException extends RuntimeException
    {
        private final Node node;

        public InvalidLayoutException(String msg, Node node)
        {
            super(msg);
            this.node = node;
        }

        /** @return the invalid Node. */
        public Node getNode()
        { return node; }
    }

    private void throwInvalidLayout(String msg, Node node)
    {
        throw new InvalidLayoutException(msg, node);
    }

    private void checkLayout(Node root)
    {
        if (root.isSplit())
        {
            Leaf split = (Leaf) root;
            if (split.getChildren().size() <= 2)
            {
                // jbw I don't care throwInvalidLayout("Split must have > 2 children", root);
            }
            Iterator<Node> splitChildren = split.getChildren().iterator();
            double weight = 0.0;
            while (splitChildren.hasNext())
            {
                Node splitChild = splitChildren.next();
                if (splitChild instanceof Divider)
                {
                    throwInvalidLayout("expected a Split or Leaf Node", splitChild);
                }
                if (splitChildren.hasNext())
                {
                    Node dividerChild = splitChildren.next();
                    if (!(dividerChild instanceof Divider))
                    {
                        throwInvalidLayout("expected a Divider Node", dividerChild);
                    }
                }
                weight += splitChild.getWeight();

                checkLayout(splitChild);
            }
            if (weight > 1.0001)
            {
                throwInvalidLayout("Split children's total weight > 1.0", root);
            }

        }
    }

    /**
     * Compute the bounds of all of the Split/Divider/Leaf Nodes in
     * the layout model, and then set the bounds of each child component
     * with a matching Leaf Node.
     */
    public void layoutContainer(Container parent)
    {
        checkLayout(getRoot());
        Insets insets = parent.getInsets();
        Dimension size = parent.getSize();
        int width = size.width - (insets.left + insets.right);
        int height = size.height - (insets.top + insets.bottom);
        Rectangle bounds = new Rectangle(insets.left, insets.top, width, height);
        layout1(getRoot(), bounds);
        layout2(getRoot(), bounds);
    }


    private Divider dividerAt(Node root, int x, int y)
    {
        if (root.type == DIVIDER)
        {
            Divider divider = (Divider) root;
            return (divider.getBounds().contains(x, y)) ? divider : null;
        }
        else if (root.isSplit())
        {
            Leaf split = (Leaf) root;
            for (Node child : split.getChildren())
            {
                if (child.getBounds().contains(x, y))
                {
                    return dividerAt(child, x, y);
                }
            }
        }
        return null;
    }

    /**
     * Return the Divider whose bounds contain the specified
     * point, or null if there isn't one.
     *
     * @param x x coordinate
     * @param y y coordinate
     * @return the Divider at x,y
     */
    public Divider dividerAt(int x, int y)
    {
        return dividerAt(getRoot(), x, y);
    }

    private boolean nodeOverlapsRectangle(Node node, Rectangle r2)
    {
        Rectangle r1 = node.getBounds();
        return
            (r1.x <= (r2.x + r2.width)) && ((r1.x + r1.width) >= r2.x) &&
                (r1.y <= (r2.y + r2.height)) && ((r1.y + r1.height) >= r2.y);
    }

    private java.util.List<Divider> dividersThatOverlap(Node root, Rectangle r)
    {
        if (nodeOverlapsRectangle(root, r) && (root.isSplit()))
        {
            java.util.List<Divider> dividers = new ArrayList();
            for (Node child : ((Leaf) root).getChildren())
            {
                if (child.type == DIVIDER)
                {
                    if (nodeOverlapsRectangle(child, r))
                    {
                        dividers.add((Divider) child);
                    }
                }
                else if (child.isSplit())
                {
                    dividers.addAll(dividersThatOverlap(child, r));
                }
            }
            return dividers;
        }
        else
        {
            return Collections.emptyList();
        }
    }

    /**
     * Return the Dividers whose bounds overlap the specified
     * Rectangle.
     *
     * @param r target Rectangle
     * @return the Dividers that overlap r
     * @throws IllegalArgumentException if the Rectangle is null
     */
    public java.util.List<Divider> dividersThatOverlap(Rectangle r)
    {
        if (r == null)
        {
            throw new IllegalArgumentException("null Rectangle");
        }
        return dividersThatOverlap(getRoot(), r);
    }


    /** Base class for the nodes that model a MultiSplitLayout. */
    public static abstract class Node
    {
        protected Leaf parent = null;
        protected String name = "";
        private Rectangle bounds = new Rectangle();
        private double weight = 0.0;
        private Component component;
        protected byte type = NONE;

        /**
         * Returns the Split parent of this Node, or null.
         *
         * @return the value of the parent property.
         * @see #setParent
         */
        public Leaf getParent()
        { return parent; }

        /**
         * Set the value of this Node's parent property.  The default
         * value of this property is null.
         *
         * @param parent a Split or null
         * @see #getParent
         */
        public void setParent(Leaf parent)
        {
            this.parent = parent;
        }

        public boolean isSplit()
        {
            return type == ROW || type == COL;
        }

        public void delete()
        {
            parent.removeChild(this);
        }

        /**
         * Returns the bounding Rectangle for this Node.
         *
         * @return the value of the bounds property.
         * @see #setBounds
         */
        public Rectangle getBounds()
        {
            return new Rectangle(this.bounds);
        }

        /**
         * Set the bounding Rectangle for this node.  The value of
         * bounds may not be null.  The default value of bounds
         * is equal to <code>new Rectangle(0,0,0,0)</code>.
         *
         * @param bounds the new value of the bounds property
         * @throws IllegalArgumentException if bounds is null
         * @see #getBounds
         */
        public void setBounds(Rectangle bounds)
        {
            if (bounds == null)
            {
                throw new IllegalArgumentException("null bounds");
            }
            this.bounds = new Rectangle(bounds);
        }

        /**
         * Value between 0.0 and 1.0 used to compute how much space
         * to add to this sibling when the layout grows or how
         * much to reduce when the layout shrinks.
         *
         * @return the value of the weight property
         * @see #setWeight
         */
        public double getWeight()
        { return weight; }

        /**
         * The weight property is a between 0.0 and 1.0 used to
         * compute how much space to add to this sibling when the
         * layout grows or how much to reduce when the layout shrinks.
         * If rowLayout is true then this node's width grows
         * or shrinks by (extraSpace * weight).  If rowLayout is false,
         * then the node's height is changed.  The default value
         * of weight is 0.0.
         *
         * @param weight a double between 0.0 and 1.0
         * @throws IllegalArgumentException if weight is not between 0.0 and 1.0
         * @see #getWeight
         */
        public void setWeight(double weight)
        {
            if ((weight < 0.0) || (weight > 1.0))
            {
                throw new IllegalArgumentException("invalid weight");
            }
            this.weight = weight;
        }

        private Node siblingAtOffset(int offset)
        {
            Leaf parent = getParent();
            if (parent == null) { return null; }
            java.util.List<Node> siblings = parent.getChildren();
            int index = siblings.indexOf(this);
            if (index == -1) { return null; }
            index += offset;
            return ((index > -1) && (index < siblings.size())) ? siblings.get(index) : null;
        }

        /**
         * Return the Node that comes after this one in the parent's
         * list of children, or null.  If this node's parent is null,
         * or if it's the last child, then return null.
         *
         * @return the Node that comes after this one in the parent's list of children.
         * @see #previousSibling
         * @see #getParent
         */
        public Node nextSibling()
        {
            return siblingAtOffset(+1);
        }

        /**
         * Return the Node that comes before this one in the parent's
         * list of children, or null.  If this node's parent is null,
         * or if it's the last child, then return null.
         *
         * @return the Node that comes before this one in the parent's list of children.
         * @see #nextSibling
         * @see #getParent
         */
        public Node previousSibling()
        {
            return siblingAtOffset(-1);
        }

        /**
         * Return the Leaf's name.
         *
         * @return the value of the name property.
         * @see #setName
         */
        public String getName()
        { return name; }

        /**
         * Set the value of the name property.  Name may not be null.
         *
         * @param name value of the name property
         * @throws IllegalArgumentException if name is null
         */
        public void setName(String name)
        {
            if (name == null)
            {
                throw new IllegalArgumentException("name is null");
            }
            this.name = name;
        }

        public Component getComponent()
        {
            return component;
        }

        public void setComponent(Component component)
        {
            this.component = component;
        }

        public void removeComponent()
        {
            component = null;
        }

        public void replaceComponent(Component component)
        {
            this.component = component;
        }
    }

    /**
     * Defines a vertical or horizontal subdivision into two or more
     * tiles.
     */
    public class Leaf extends Node
    {
        ArrayList<Node> children = new ArrayList();

        public Leaf()
        {
            name = "root";
            type = ROOT;
        }

        public Leaf(byte type)
        {
            name = "root";
            this.type = type;
        }

        public Leaf(String name)
        {
            this.name = name;
            this.type = LEAF;
        }

        public Leaf(String name, byte type)
        {
            this.name = name;
            this.type = type;
        }

        public void setLeafType(byte type)
        {
            if(type != ROW && type != COL)
            {
                StsException.systemError(this, "setLeafType", "Attempting to set splitType to illegal vaue: " + type);
                return;
            }
            else if(this.type == type)
                return;
            else if(isSplit())
                StsException.systemError(this, "setLeafType", "Leaf is already split as: " + this.type + ". Can't split as: " + type);
            else
                this.type = type;
        }

        /**
         * Returns true if the this Split's children are to be
         * laid out in a row: all the same height, left edge
         * equal to the previous Node's right edge.  If false,
         * children are laid on in a column.
         *
         * @return the value of the rowLayout property.
         */
        public boolean isRowLayout()
        {
            return type == ROW;
        }

        /**
         * Returns this Split node's children.  The returned value
         * is not a reference to the Split's internal list of children
         *
         * @return the value of the children property.
         */
        public ArrayList<Node> getChildren()
        {
            return new ArrayList<Node>(children);
        }

        /**
         * Set's the children property of this Leaf node.  The parent
         * of each new child is set to this Leaf node, and the parent
         * of each old child (if any) is set to null.  This method
         * defensively copies the incoming List.  Default value is
         * an empty List.
         *
         * @throws IllegalArgumentException if children is null
         * @see #getChildren
         */
        /*
        public void setChildren(java.util.List<Node> children)
        {
            if (children == null)
            {
                throw new IllegalArgumentException("children must be a non-null List");
            }
            for (Node child : this.children)
            {
                child.setParent(null);
            }
            this.children = new ArrayList<Node>(children);
            for (Node child : this.children)
            {
                child.setParent(this);
            }
        }
       */
        public ArrayList addLeafChildren(String[] leafNames)
        {
            ArrayList leaves = new ArrayList();
            for (String leafName : leafNames)
            {
                Leaf child = new Leaf(leafName);
                leaves.add(child);
                addLeafChild(child);
            }
            return leaves;
        }

        public Leaf addLeafChild(String leafName)
        {
            Leaf child = new Leaf(leafName);
            addLeafChild(child);
            return child;
        }

        private void addLeafChild(Leaf child)
        {
            if (children.size() == 0)
            {
                child.setParent(this);
                children.add(child);
            }
            else
            {
                Divider divider = new Divider();
                divider.setParent(this);
                children.add(divider);
                child.setParent(this);
                children.add(child);
				child.setWeight(.05);
            }
            leafList.add(child);
//            setChildren(children);
        }
        /**
         * Convenience method that returns the last child whose weight
         * is > 0.0.
         *
         * @return the last child whose weight is > 0.0.
         * @see #getChildren
         * @see Node#getWeight
         */
        public final Node lastWeightedChild()
        {
            java.util.List<Node> children = getChildren();
            Node weightedChild = null;
            for (Node child : children)
            {
                if (child.getWeight() > 0.0)
                {
                    weightedChild = child;
                }
            }
            return weightedChild;
        }

        public void replaceChild(Node oldChild, Node newChild)
        {
            boolean replaced = false;
            ArrayList<Node> newChildren = new ArrayList();
            for (Node child : children)
            {
                if (child == oldChild)
                {
                    newChildren.add(newChild);
                    replaced = true;
                }
                else
                    newChildren.add(child);
            }
            children = newChildren;
            if (!replaced) System.err.println("Failed to replace child in list");
        }

        public void removeChild(Node removeChild)
        {
            if(children.size() == 1)
            {
                StsException.systemError(this, "removeChild", "nChildren is one. Should never happen.");
                return;
            }

//            removeChild.removeComponent();
            int childIndex = children.indexOf(removeChild);
            if(childIndex == -1)
            {
                StsException.systemError(this, "removeChild", "Failed to find child: " + removeChild.getName() +
                " in children list of split node: " + getName());
                return;
            }

            if(childIndex > 0)
            {
                children.remove(childIndex);
                children.remove(childIndex-1); // remove previous divider
            }
            else
            {
                children.remove(childIndex+1); // remove next divider
                children.remove(childIndex);

            }
            int nChildren = children.size();
            if(nChildren == 1)
            {
                Leaf grandParent = getParent();
				if (grandParent != null) // jbw
				{
					grandParent.replaceChild(this, children.get(0));
					children.remove(0);
				}
            }
        }

        public String toString()
        {
            int nChildren = getChildren().size();
            StringBuffer sb = new StringBuffer("MultiSplitLayout.Split");
            sb.append(isRowLayout() ? " ROW [" : " COLUMN [");
            sb.append(nChildren + ((nChildren == 1) ? " child" : " children"));
            sb.append("] ");
            sb.append(getBounds());
            return sb.toString();
        }
    }

    /** Models a single vertical/horiztonal divider. */
    public class Divider extends Node
    {

        public Divider()
        {
            type = DIVIDER;
        }
        /**
         * Convenience method, returns true if the Divider's parent
         * is a Split row (a Split with isRowLayout() true), false
         * otherwise. In other words if this Divider's major axis
         * is vertical, return true.
         *
         * @return true if this Divider is part of a Split row.
         */
        public final boolean isVertical()
        {
            Leaf parent = getParent();
            return (parent != null) ? parent.isRowLayout() : false;
        }

        /**
         * Dividers can't have a weight, they don't grow or shrink.
         *
         * @throws UnsupportedOperationException
         */
        public void setWeight(double weight)
        {
            throw new UnsupportedOperationException();
        }

        public String toString()
        {
            return "MultiSplitLayout.Divider " + getBounds().toString();
        }
    }

    public void addRootLeaf(String leafName)
    {
        root = new Leaf(LEAF);
        root.setName(leafName);
		root.setWeight(.05);
        leafList.add(root);
    }

    public void addRootRowSplit(String[] leafNames)
    {
        Leaf leaf = new Leaf(ROW);
        root = leaf;
        leafList.add(leaf);
        ArrayList leaves = leaf.addLeafChildren(leafNames);
//        leafList.addAll(leaves);
    }

    public void addEqualRootRowSplit(String[] leafNames)
    {
        addRootRowSplit(leafNames);
        Leaf leaf;
        int nLeafs = leafNames.length;
        double weight = 1.0/nLeafs;
        for(int n = 0; n < leafNames.length; n++)
        {
            leaf = getLeaf(leafNames[n]);
            leaf.setWeight(weight);
        }
    }

    public Leaf addRootRowSplit(String leafName)
    {
        if(!root.isSplit())
        {
            Leaf rootLeaf = root;
            Leaf newRootLeaf = new Leaf(ROW);
            root = newRootLeaf;
            leafList.add(newRootLeaf);
            newRootLeaf.addLeafChild(rootLeaf);
            Leaf newLeaf = new Leaf(leafName);
            newRootLeaf.addLeafChild(newLeaf);
            return newLeaf;
        }
        else
            return root.addLeafChild(leafName);
    }

    public void addRowSplit(String parentName, String[] leafNames)
    {
        addSplit(parentName, ROW, leafNames);
    }

    public void addRootColSplit(String[] leafNames)
    {
        Leaf leaf = new Leaf(COL);
        root = leaf;
        leafList.add(leaf);
        ArrayList leaves = leaf.addLeafChildren(leafNames);
//        leafList.addAll(leaves);
    }

    public void addColSplit(String parentName, String[] leafNames)
    {
        addSplit(parentName, COL, leafNames);
    }
	public void addEqualRootColSplit(String[] leafNames)
	{
		addRootColSplit(leafNames);
		Leaf leaf;
		int nLeafs = leafNames.length;
		double weight = 1.0/nLeafs;
		for(int n = 0; n < leafNames.length; n++)
		{
			leaf = getLeaf(leafNames[n]);
			leaf.setWeight(weight);
		}
    }
    public void addSplit(String parentName, byte type, String[] leafNames)
    {
        for(Leaf split : leafList)
        {
            if (split.getName().equals(parentName))
            {
                split.setLeafType(type);
                ArrayList leaves = split.addLeafChildren(leafNames);
//                leafList.addAll(leaves);
                return;
            }
        }
    }

    public void addChild(String parentName, String leafName)
    {
       for(Leaf parent : leafList)
        {
            if (parent.getName().equals(parentName))
            {
                ArrayList children = parent.getChildren();
                Leaf leaf = new Leaf(leafName);
                children.add(leaf);
                leafList.add(leaf);
                return;
            }
        }
        StsException.systemError(this, "addChild", "Failed to find parent: " + parentName + " in order to add child " + leafName);
    }

    public void addRootChild(String leafName)
    {
        root.addLeafChild(leafName);
    }

    public void setWeight(String name, double weight)
    {
        Leaf leaf = getLeaf(name);
        if(leaf == null)
        {
            StsException.systemError(this, "setWeight", "Failed to find leaf: " + name);
            return;
        }
        leaf.setWeight(weight);
    }

    public Leaf getLeaf(String name)
    {
        for(Leaf leaf : leafList)
            if(leaf.getName().equals(name))
                return leaf;
        return null;
    }

    private void constructLeafList()
    {
        LeafIterator iterator = new LeafIterator();
        leafList = iterator.list;
    }

    class LeafIterator implements Iterator
    {
        ArrayList<Leaf> list;
        private Iterator iterator;

        LeafIterator()
        {
            list = LeafIterator(root);
        }

        private ArrayList LeafIterator(Node node)
        {
            list = new ArrayList();
            if (node.type != DIVIDER)
            {
                Leaf leaf = (Leaf)node;
                list.add(leaf);
                ArrayList<Node> children = leaf.getChildren();
                iterator = children.iterator();
                while (iterator.hasNext())
                    list.addAll(LeafIterator((Leaf)iterator.next()));
            }

            return list;
        }

        public boolean hasNext()
        {
            return iterator.hasNext();
        }

        public synchronized Object next()
        {
            return iterator.next();
        }

        public void remove()
        {
            throw new UnsupportedOperationException();
        }
    }

    static public void main(String[] args)
    {
        StsMultiSplitPaneLayout layout = new StsMultiSplitPaneLayout();
        layout.addRootRowSplit(new String[]{"left", "center", "right"});
    }
}
