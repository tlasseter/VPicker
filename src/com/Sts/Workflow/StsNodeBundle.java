package com.Sts.Workflow;

import com.Sts.Utilities.*;
/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

// A prefix of nodes with an associated constraint used to determine if input requirements for
// a node are satisfied.

public class StsNodeBundle
{
    StsWorkflowTreeNode[] nodes;
    int constraint;

    static public final int NONE_REQUIRED = 0;
    static public final int ONE_REQUIRED = 1;
    static public final int ALL_REQUIRED = 2;

    public StsNodeBundle(StsWorkflowTreeNode[] nodes, int constraint)
    {
        this.nodes = nodes;
        if(nodes == null) nodes = new StsWorkflowTreeNode[0];
        this.constraint = constraint;
    }

    public StsNodeBundle(StsWorkflowTreeNode node)
    {
        if(node == null) nodes = new StsWorkflowTreeNode[0];
        else             nodes = new StsWorkflowTreeNode[] { node };
        constraint = ALL_REQUIRED;
    }

    public void addNodetoBundle(StsWorkflowTreeNode node)
    {
        if(nodes == null) nodes = new StsWorkflowTreeNode[0];
        else
            nodes = (StsWorkflowTreeNode[])StsMath.arrayAddArray(nodes, node);
    }

    public void addOutputNode(StsWorkflowTreeNode outNode)
    {
        int nNodes = nodes.length;
        for(int n = 0; n < nNodes; n++)
            nodes[n].addOutputNode(outNode);
    }

    protected boolean constraintSatisfied()
    {
        StsWorkflowTreeNode node;

        if(constraint == NONE_REQUIRED) return true;

        boolean satisfied = true;
        int nNodes = nodes.length;
        int nSatisfied = 0;
        for(int n = 0; n < nNodes; n++)
            if(nodes[n].getActionStatus() > StsWorkflowTreeNode.CAN_START) nSatisfied++;

        if(constraint == ONE_REQUIRED) return nSatisfied > 0;
        else return nSatisfied == nNodes;
    }

    protected boolean adjustedState(String actionClassName, int actionStatus)
    {
        int nNodes = nodes.length;
        for(int n = 0; n < nNodes; n++)
            if(nodes[n].actionClassName == actionClassName)
            {
                nodes[n].setActionStatus(actionStatus);
                return true;
            }
        return false;
    }

    protected void setActionStatus(int actionStatus)
    {
        for(int n = 0; n < nodes.length; n++)
            nodes[n].setActionStatus(actionStatus);
    }
}
