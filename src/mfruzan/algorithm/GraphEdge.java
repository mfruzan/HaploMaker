/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mfruzan.algorithm;

/**
 *
 * @author a1195806
 */
public class GraphEdge {
     public TreeNode from = null;
    public TreeNode to = null;
    public boolean isDummy = false;
    public GraphEdge(TreeNode f, TreeNode t){
        from = f;
        to =t;
    }
}
