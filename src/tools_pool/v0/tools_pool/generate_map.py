import sys
import os
import json
import random
from typing import List, Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, Polygon
from matplotlib.text import Text
import matplotlib.patheffects as path_effects

# Add the parent directory to the sys.path to find common_utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from common_utils.sequence import SequenceRecord, Feature, load_sequence_from_json

class PlasmidMapper:
    """
    根据 SequenceRecord 数据生成专业级质粒图谱的类。
    """
    def __init__(self, record: SequenceRecord):
        """
        初始化质粒图谱生成器。
        :param record: 一个符合 Schema 定义的 SequenceRecord 对象。
        """
        if not record.get("circular", False):
            raise ValueError("该可视化工具目前仅支持环状质粒。")

        self.record = record
        self.sequence_length = record.get("length", 0)
        self.title = record.get("name", "Unnamed Plasmid")
        self.features_to_plot = []
        self.sites_to_plot = []
        self.scale_marks = []
        self._populate_from_record()

    def _populate_from_record(self):
        """从 SequenceRecord 中提取并转换绘图所需的数据。"""
        features = self.record.get("features", [])
        
        # 1. 分类特征：酶切位点 vs. 其他特征
        for feat in features:
            feat_type = feat.get("type", "misc_feature")
            label = feat.get("qualifiers", {}).get("label", feat.get("qualifiers", {}).get("name", ""))
            start = feat.get("start")
            end = feat.get("end")

            if not all([label, start, end]):
                continue # 跳过信息不完整的特征

            if feat_type == "restriction_site":
                self.sites_to_plot.append({
                    "position": start,
                    "name": label,
                    "color": "#808080", # 默认为灰色
                    "size": 0.1,
                    "label_position": True
                })
            else:
                strand = feat.get("strand", 1)
                direction = "clockwise" if strand == 1 else "counter-clockwise"
                self.features_to_plot.append({
                    "start": start,
                    "end": end,
                    "label": label,
                    "color": "#66ccff", # 默认颜色，后续会根据类型自动匹配
                    "height": 0.05,
                    "plot_label": True,
                    "feature_type": feat_type,
                    "arrow": True, # 默认显示箭头
                    "direction": direction
                })
        
        # 2. 自动生成刻度标记
        if self.sequence_length > 0:
            step = 1000 if self.sequence_length >= 5000 else 500
            for i in range(0, self.sequence_length, step):
                self.scale_marks.append({"position": i, "label": str(i)})

    def generate_map(self, output_path="plasmid_map.png", dpi=300, figsize=(10, 10)):
        """
        生成并保存专业级环状质粒图谱。
        """
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_aspect('equal')
        ax.axis('off')
        
        COLOR_SCHEME = {
            "promoter": "#FFD700", "terminator": "#FF6347", "resistance": "#32CD32",
            "origin": "#1E90FF", "gene": "#87CEEB", "cds": "#87CEEB", "tag": "#BA55D3",
            "enhancer": "#FFA07A", "misc_feature": "#A9A9A9", "default": "#66ccff"
        }
        
        center = (0.5, 0.5)
        outer_radius = 0.45
        inner_radius = outer_radius - 0.01
        
        # 绘制质粒骨架
        ax.add_patch(Circle(center, outer_radius, fill=False, edgecolor='#0f0700', linewidth=2))
        ax.add_patch(Circle(center, inner_radius, fill=False, edgecolor='#0f0700', linewidth=2))
        
        # 绘制刻度
        for mark in self.scale_marks:
            angle = np.radians(360 - (360 * mark["position"] / self.sequence_length) + 90)
            ax.plot(
                [center[0] + inner_radius * np.cos(angle), center[0] + (inner_radius - 0.015) * np.cos(angle)],
                [center[1] + inner_radius * np.sin(angle), center[1] + (inner_radius - 0.015) * np.sin(angle)],
                'k-', linewidth=1
            )
            if mark["label"] == "0": # 仅在0点标记数字
                 ax.text(
                    center[0] + (inner_radius - 0.03) * np.cos(angle),
                    center[1] + (inner_radius - 0.03) * np.sin(angle),
                    mark["label"], ha='center', va='center', fontsize=10, fontweight='bold'
                )

        # 绘制特征
        for feat in self.features_to_plot:
            start, end = feat["start"], feat["end"]
            start_angle_deg = 360 - (360 * start / self.sequence_length) + 90
            end_angle_deg = 360 - (360 * end / self.sequence_length) + 90

            if start > end: # 处理跨越原点的特征
                end_angle_deg += 360

            color = COLOR_SCHEME.get(feat["feature_type"], COLOR_SCHEME["default"])
            arc_radius = inner_radius - 0.04
            
            arc = Arc(
                center, 2 * arc_radius, 2 * arc_radius,
                theta1=end_angle_deg, theta2=start_angle_deg,
                color=color, linewidth=10, alpha=0.8
            )
            ax.add_patch(arc)

            # 绘制特征标签 (简化版，放置在弧线中点外部)
            if feat["plot_label"]:
                mid_angle_rad = np.radians((start_angle_deg + end_angle_deg) / 2)
                label_radius = arc_radius + 0.05
                text_x = center[0] + label_radius * np.cos(mid_angle_rad)
                text_y = center[1] + label_radius * np.sin(mid_angle_rad)
                ax.text(text_x, text_y, feat["label"], ha='center', va='center', fontsize=9, fontweight='bold')

            # 绘制箭头
            if feat.get("arrow", False):
                arrow_angle_deg = end_angle_deg if feat["direction"] == "clockwise" else start_angle_deg
                arrow_rad = np.radians(arrow_angle_deg)
                
                arrow_base_x = center[0] + arc_radius * np.cos(arrow_rad)
                arrow_base_y = center[1] + arc_radius * np.sin(arrow_rad)
                
                arrow_size = 0.015
                dx = np.sin(arrow_rad) * arrow_size
                dy = -np.cos(arrow_rad) * arrow_size
                
                if feat["direction"] == "counter-clockwise":
                    dx, dy = -dx, -dy

                arrow_points = [
                    (arrow_base_x + dx, arrow_base_y + dy),
                    (arrow_base_x - 0.5*dx + 0.8*dy, arrow_base_y - 0.5*dy - 0.8*dx),
                    (arrow_base_x - 0.5*dx - 0.8*dy, arrow_base_y - 0.5*dy + 0.8*dx),
                ]
                ax.add_patch(Polygon(arrow_points, closed=True, color=color, zorder=10))

        # 绘制酶切位点
        for site in self.sites_to_plot:
            angle_rad = np.radians(360 - (360 * site["position"] / self.sequence_length) + 90)
            site_x = center[0] + outer_radius * np.cos(angle_rad)
            site_y = center[1] + outer_radius * np.sin(angle_rad)
            
            label_radius = outer_radius + 0.05
            label_x = center[0] + label_radius * np.cos(angle_rad)
            label_y = center[1] + label_radius * np.sin(angle_rad)

            ax.plot([site_x, label_x], [site_y, label_y], color=site["color"], linestyle='-', linewidth=0.8)
            label_text = f"{site['name']} ({site['position']})"
            ax.text(label_x, label_y, label_text, ha='center', va='center', fontsize=9, color=site["color"])

        # 添加标题和尺寸
        ax.set_title(self.title, fontsize=16, fontweight='bold', y=1.0)
        ax.text(0.5, 0.5, f"{self.sequence_length} bp", ha='center', va='center', fontsize=12)
        
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close()

def generate_map(json_path: str, output_path: Optional[str] = None) -> Dict[str, str]:
    """
    从 SequenceRecord JSON 文件生成质粒图谱可视化。

    Args:
        json_path (str): 输入的 SequenceRecord JSON 文件路径。
        output_path (Optional[str]): 输出图像文件的路径。如果未提供，
                                     将根据输入文件名在同一目录下生成一个 .png 文件。

    Returns:
        Dict[str, str]: 包含输出文件路径的字典。
    """
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"输入文件未找到: {json_path}")

    if output_path is None:
        base_name = os.path.splitext(os.path.basename(json_path))[0]
        output_path = os.path.join(os.path.dirname(json_path), f"{base_name}_map.png")

    record = load_sequence_from_json(json_path)
    
    mapper = PlasmidMapper(record)
    mapper.generate_map(output_path=output_path)
    
    abs_output_path = os.path.abspath(output_path)
    print(f"质粒图谱已生成: {abs_output_path}")
    
    return {"output_path": abs_output_path}

if __name__ == "__main__":
    # 确保示例 JSON 文件存在
    example_json_path = os.path.join("data", "temp", "pcDNA3.1(-).json")
    
    if not os.path.exists(example_json_path):
        print(f"错误: 示例文件 '{example_json_path}' 不存在。")
        print("请先运行 `common_utils/loadSequence.py` 来生成该文件。")
    else:
        print(f"正在使用示例文件 '{example_json_path}' 生成图谱...")
        result = generate_map(example_json_path)
        print(f"图谱保存在: {result['output_path']}")
