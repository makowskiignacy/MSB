#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lung Nodules Classification with XAI - Mac M-series Optimized Version
Optimized for local execution on Apple Silicon with file-based result saving
"""

import sys
import os
import json
import pickle
from datetime import datetime
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import torch
from torch.utils.data import DataLoader
from torchvision.transforms import v2
from torch import nn
import numpy as np
import matplotlib.pyplot as plt
import gdown
from transformers import ViTConfig, ViTModel
from medmnist import NoduleMNIST3D, PathMNIST
from medmnist.info import INFO as MEDMNIST_INFO

# Configure for Apple Silicon
if torch.backends.mps.is_available():
    DEVICE = "mps"  # Apple Silicon GPU
    print("Using Apple Silicon GPU (MPS)")
elif torch.cuda.is_available():
    DEVICE = "cuda"
    print("Using CUDA GPU")
else:
    DEVICE = "cpu"
    print("Using CPU")

# Set up results directory
RESULTS_DIR = Path("./results")
RESULTS_DIR.mkdir(exist_ok=True)

# Create subdirectories for different types of results
(RESULTS_DIR / "models").mkdir(exist_ok=True)
(RESULTS_DIR / "visualizations").mkdir(exist_ok=True)
(RESULTS_DIR / "predictions").mkdir(exist_ok=True)
(RESULTS_DIR / "explanations").mkdir(exist_ok=True)
(RESULTS_DIR / "logs").mkdir(exist_ok=True)

print(f"Results will be saved to: {RESULTS_DIR.absolute()}")

def log_message(message, log_file="execution_log.txt"):
    """Log messages to file with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_path = RESULTS_DIR / "logs" / log_file
    with open(log_path, "a") as f:
        f.write(f"[{timestamp}] {message}\n")
    print(f"[{timestamp}] {message}")

def save_json(data, filename, subdir=""):
    """Save data as JSON file"""
    if subdir:
        save_path = RESULTS_DIR / subdir / filename
    else:
        save_path = RESULTS_DIR / filename
    
    with open(save_path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
    log_message(f"Saved JSON to {save_path}")

def save_numpy(data, filename, subdir=""):
    """Save numpy array"""
    if subdir:
        save_path = RESULTS_DIR / subdir / filename
    else:
        save_path = RESULTS_DIR / filename
    
    np.save(save_path, data)
    log_message(f"Saved numpy array to {save_path}")

def save_figure(fig, filename, subdir="visualizations", dpi=300):
    """Save matplotlib figure"""
    save_path = RESULTS_DIR / subdir / filename
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    log_message(f"Saved figure to {save_path}")

def download_weights(url, output_dir, filename):
    """Downloads weights with progress tracking"""
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)

    if not os.path.exists(output_path):
        log_message(f"Downloading weights to {output_path}...")
        gdown.download(url, output_path)
        log_message("Download completed")
    else:
        log_message(f"Weights already exist at {output_path}. Skipping download.")
    
    return output_path

class DINO(nn.Module):
    """DINO Transformer model optimized for Apple Silicon"""
    
    def __init__(self):
        super().__init__()
        # Backbone with eager attention for compatibility
        config = ViTConfig.from_pretrained('facebook/dino-vits8', attn_implementation="eager")
        self.backbone = ViTModel(config)
        # Classification head
        self.head = torch.nn.Linear(384, 1)

    def forward(self, x: torch.Tensor, output_attentions: bool = False):
        out = self.backbone(x, output_attentions=output_attentions)
        x = out["pooler_output"]
        x = self.head(x)
        if output_attentions:
            att = out["attentions"]
            return x, att
        else:
            return x

def setup_lung_nodule_data():
    """Setup NoduleMNIST data with optimizations for Mac"""
    log_message("Setting up NoduleMNIST dataset...")
    
    def take_middle_slice(inpt: np.ndarray):
        """Extract middle slice and convert to 3-channel format"""
        inpt = inpt.squeeze()
        X, Y, Z = inpt.shape
        slice_ = inpt[:, :, Z//2]
        slice_ = torch.Tensor(slice_).unsqueeze(dim=0).repeat(3, 1, 1)
        return slice_

    # Optimized transforms for Apple Silicon
    transforms = v2.Compose([
        v2.Lambda(take_middle_slice),
        v2.Resize(size=(224, 224), antialias=True)
    ])

    normalize = v2.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])

    # Dataset setup
    data_dir = "./example_data"
    os.makedirs(data_dir, exist_ok=True)

    # Smaller batch sizes for stability on Mac
    batch_size = 16 if DEVICE == "mps" else 32

    ref_set = NoduleMNIST3D(root=data_dir, split="val", size=64, transform=transforms, download=True)
    inf_set = NoduleMNIST3D(root=data_dir, split="test", size=64, transform=transforms)

    ref_loader = DataLoader(ref_set, batch_size=batch_size, shuffle=False, 
                           num_workers=0 if DEVICE == "mps" else 2, pin_memory=False)
    inf_loader = DataLoader(inf_set, batch_size=6, shuffle=True, 
                           num_workers=0 if DEVICE == "mps" else 2, pin_memory=False)

    class_names = ["benign", "malignant"]
    logit2name = {0: "benign", 1: "malignant"}

    # Save dataset info
    dataset_info = {
        "reference_samples": len(ref_set),
        "inference_samples": len(inf_set),
        "batch_size": batch_size,
        "class_names": class_names,
        "device": DEVICE
    }
    save_json(dataset_info, "dataset_info.json")

    log_message(f"Dataset setup complete. Ref: {len(ref_set)}, Test: {len(inf_set)}")
    
    return ref_loader, inf_loader, normalize, class_names, logit2name

def load_pretrained_model():
    """Load and setup the pretrained DINO model"""
    log_message("Loading pretrained DINO model...")
    
    # Download weights
    url = "https://drive.google.com/uc?id=1xUevCbvII5yXDxVxb7bR65CPmgz2sGQA"
    weights_path = download_weights(url, "tuned_models", "lidc_dino_s8.pth")
    
    # Initialize and load model
    model = DINO()
    
    # Load with proper device mapping for Apple Silicon
    if DEVICE == "mps":
        # MPS doesn't support all operations, load on CPU first
        state_dict = torch.load(weights_path, map_location="cpu", weights_only=True)
        model.load_state_dict(state_dict)
        model = model.to(DEVICE)
    else:
        state_dict = torch.load(weights_path, map_location=DEVICE, weights_only=True)
        model.load_state_dict(state_dict)
        model = model.to(DEVICE)
    
    model.eval()
    
    log_message("Model loaded successfully")
    return model

def visualize_samples_and_save(samples, labels, class_names, filename_prefix="samples"):
    """Visualize samples and save to file"""
    log_message(f"Creating visualization: {filename_prefix}")
    
    fig, axes = plt.subplots(1, 5, figsize=(20, 5))
    for i in range(min(5, len(samples))):
        image = samples[i].permute(1, 2, 0).cpu().numpy()
        axes[i].imshow(image, cmap='gray')
        axes[i].set_title(f"Label: {class_names[labels[i].item()]}")
        axes[i].axis('off')

    plt.tight_layout()
    save_figure(fig, f"{filename_prefix}.png")
    
    return fig

def make_predictions_and_save(model, samples, labels, class_names, logit2name):
    """Make predictions and save results"""
    log_message("Making predictions...")
    
    # Move to device with proper handling for MPS
    if DEVICE == "mps":
        # Process in smaller batches for MPS stability
        samples = samples[:5].to(DEVICE)
        labels = labels[:5].to(DEVICE)
    else:
        samples = samples[:5].to(DEVICE)
        labels = labels[:5].to(DEVICE)

    # Make predictions
    with torch.no_grad():
        logits = model(samples)
        predictions = torch.sigmoid(logits).round().squeeze().cpu().numpy()
        probabilities = torch.sigmoid(logits).squeeze().cpu().numpy()

    # Convert to interpretable format
    if predictions.ndim == 0:
        predictions = [predictions.item()]
        probabilities = [probabilities.item()]
    
    predicted_classes = [logit2name[int(pred)] for pred in predictions]
    true_classes = [class_names[labels[i].item()] for i in range(len(labels))]

    # Save predictions
    results = {
        "predictions": {
            "predicted_classes": predicted_classes,
            "true_classes": true_classes,
            "probabilities": probabilities.tolist() if hasattr(probabilities, 'tolist') else [probabilities],
            "accuracy": sum(p == t for p, t in zip(predicted_classes, true_classes)) / len(predicted_classes)
        },
        "model_info": {
            "device": DEVICE,
            "model_type": "DINO_ViT",
            "task": "binary_classification"
        }
    }
    
    save_json(results, "predictions_results.json", "predictions")
    
    # Print results
    for i, pred_class in enumerate(predicted_classes):
        log_message(f"Sample {i + 1}: Predicted={pred_class}, True={true_classes[i]}, Prob={probabilities[i]:.3f}")

    return samples, predictions, probabilities, predicted_classes

def setup_xai_tools(model):
    """Setup XAI tools with error handling"""
    log_message("Setting up XAI tools...")
    
    try:
        from obzai.xai.xai_tool import CDAM, AttentionMap
        
        # CDAM tool
        cdam_tool = CDAM(
            model=model,
            mode='vanilla',
            gradient_type="from_logits",
            gradient_reduction="average",
            activation_type="sigmoid"
        )
        cdam_tool.create_hooks(layer_name="backbone.encoder.layer.11.layernorm_before")
        
        # Attention tool
        attention_tool = AttentionMap(
            model=model,
            attention_layer_id=-1,
            head=None
        )
        
        log_message("XAI tools configured successfully")
        return cdam_tool, attention_tool
        
    except ImportError as e:
        log_message(f"XAI tools not available: {e}")
        return None, None
    except Exception as e:
        log_message(f"Error setting up XAI tools: {e}")
        return None, None

def generate_attention_maps(attention_tool, samples, class_names, true_labels, pred_labels):
    """Generate and save attention maps"""
    if attention_tool is None:
        log_message("Attention tool not available, skipping...")
        return None
        
    log_message("Generating attention maps...")
    
    try:
        # Generate attention maps
        attention_maps = attention_tool.explain(samples)
        
        if attention_maps is None:
            log_message("Failed to generate attention maps")
            return None
            
        # Save raw attention maps
        save_numpy(attention_maps.cpu().numpy(), "attention_maps.npy", "explanations")
        
        # Create visualizations
        fig, axes = plt.subplots(2, 5, figsize=(20, 10))
        
        # Original images
        for i in range(5):
            original_image = samples[i].permute(1, 2, 0).cpu().numpy()
            axes[0, i].imshow(original_image, cmap='gray')
            axes[0, i].set_title(f"Sample {i + 1}")
            axes[0, i].axis('off')
        
        # Attention maps
        for i in range(5):
            attention_map = attention_maps[i].cpu().numpy()
            axes[1, i].imshow(attention_map, cmap='jet')
            axes[1, i].set_title(f"Attention {i + 1}")
            axes[1, i].axis('off')
        
        plt.tight_layout()
        save_figure(fig, "attention_maps_comparison.png", "explanations")
        
        # Overlay visualization
        fig_overlay, axes_overlay = plt.subplots(1, 5, figsize=(20, 5))
        for i in range(5):
            original_image = samples[i].permute(1, 2, 0).cpu().numpy()
            attention_map = attention_maps[i].cpu().numpy()
            
            axes_overlay[i].imshow(original_image, cmap='gray')
            axes_overlay[i].imshow(attention_map, cmap='jet', alpha=0.5)
            axes_overlay[i].set_title(f"Sample {i + 1}")
            axes_overlay[i].axis('off')
        
        plt.tight_layout()
        save_figure(fig_overlay, "attention_maps_overlay.png", "explanations")
        
        log_message("Attention maps generated and saved")
        return attention_maps
        
    except Exception as e:
        log_message(f"Error generating attention maps: {e}")
        return None

def generate_cdam_maps(cdam_tool, samples, target_indices):
    """Generate and save CDAM maps"""
    if cdam_tool is None:
        log_message("CDAM tool not available, skipping...")
        return None
        
    log_message("Generating CDAM maps...")
    
    try:
        # Generate CDAM maps
        cdam_maps = cdam_tool.explain(samples, target_idx=target_indices)
        
        if cdam_maps is None:
            log_message("Failed to generate CDAM maps")
            return None
            
        # Save raw CDAM maps
        save_numpy(cdam_maps.cpu().numpy(), "cdam_maps.npy", "explanations")
        
        # Create visualizations
        fig, axes = plt.subplots(2, 5, figsize=(20, 10))
        
        # Original images
        for i in range(5):
            original_image = samples[i].permute(1, 2, 0).cpu().numpy()
            axes[0, i].imshow(original_image, cmap='gray')
            axes[0, i].set_title(f"Sample {i + 1}")
            axes[0, i].axis('off')
        
        # CDAM maps
        cdam_abs_max = cdam_maps.abs().max()
        for i in range(5):
            cdam_map = cdam_maps[i].squeeze().cpu().numpy()
            axes[1, i].imshow(cdam_map, cmap='coolwarm', vmin=-cdam_abs_max, vmax=cdam_abs_max)
            axes[1, i].set_title(f"CDAM {i + 1}")
            axes[1, i].axis('off')
        
        plt.tight_layout()
        save_figure(fig, "cdam_maps_comparison.png", "explanations")
        
        # Overlay with histograms
        fig_detailed, axes_detailed = plt.subplots(2, 5, figsize=(20, 10), 
                                                  gridspec_kw={'height_ratios': [4, 1]})
        
        for i in range(5):
            original_image = samples[i].permute(1, 2, 0).cpu().numpy()
            cdam_map = cdam_maps[i].squeeze().cpu().numpy()
            
            # Overlay
            axes_detailed[0, i].imshow(original_image, cmap='gray')
            im = axes_detailed[0, i].imshow(cdam_map, cmap='coolwarm', alpha=0.5, 
                                          vmin=-cdam_abs_max, vmax=cdam_abs_max)
            axes_detailed[0, i].set_title(f"Sample {i + 1}")
            axes_detailed[0, i].axis('off')
            
            # Histogram
            axes_detailed[1, i].hist(cdam_map.ravel(), bins=30, color='blue', alpha=0.7)
            axes_detailed[1, i].set_title(f"Histogram {i + 1}")
            axes_detailed[1, i].set_xlabel('Value')
            axes_detailed[1, i].set_ylabel('Frequency')
        
        plt.tight_layout()
        save_figure(fig_detailed, "cdam_maps_detailed.png", "explanations")
        
        log_message("CDAM maps generated and saved")
        return cdam_maps
        
    except Exception as e:
        log_message(f"Error generating CDAM maps: {e}")
        return None

def setup_outlier_detection(ref_loader):
    """Setup outlier detection with error handling"""
    log_message("Setting up outlier detection...")
    
    try:
        from obzai.data_inspector.extractor import FirstOrderExtractor
        from obzai.data_inspector.detector import GMMDetector
        
        # Setup extractor and detector
        first_order_extrc = FirstOrderExtractor()
        gmm_detector = GMMDetector(
            extractors=[first_order_extrc], 
            n_components=3, 
            outlier_quantile=0.01, 
            show_progress=True
        )
        
        # Fit on reference data
        gmm_detector.fit(ref_loader)
        
        # Save detector
        with open(RESULTS_DIR / "models" / "outlier_detector.pkl", 'wb') as f:
            pickle.dump(gmm_detector, f)
        
        log_message("Outlier detection setup complete")
        return gmm_detector
        
    except ImportError as e:
        log_message(f"Outlier detection not available: {e}")
        return None
    except Exception as e:
        log_message(f"Error setting up outlier detection: {e}")
        return None

def run_homework_pathmnist():
    """Run the homework section for PathMNIST with file saving"""
    log_message("Starting PathMNIST homework...")
    
    try:
        # Setup PathMNIST
        dataset_name = 'pathmnist'
        info = MEDMNIST_INFO[dataset_name]
        n_classes = len(info['label'])
        class_names = [info['label'][str(i)] for i in range(n_classes)]
        
        log_message(f"PathMNIST: {n_classes} classes, Names: {class_names}")
        
        # Transforms for PathMNIST (optimized for Mac)
        transforms = v2.Compose([
            v2.ToTensor(),
            v2.Resize(size=(224, 224), antialias=True),
            v2.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
        ])
        
        # Data setup
        hw_data_dir = "./homework_data_pathmnist"
        os.makedirs(hw_data_dir, exist_ok=True)
        
        batch_size = 16 if DEVICE == "mps" else 32
        num_workers = 0 if DEVICE == "mps" else 2
        
        hw_train_dataset = PathMNIST(root=hw_data_dir, split="train", transform=transforms, download=True)
        hw_val_dataset = PathMNIST(root=hw_data_dir, split="val", transform=transforms, download=True)
        hw_test_dataset = PathMNIST(root=hw_data_dir, split="test", transform=transforms, download=True)
        
        hw_train_loader = DataLoader(hw_train_dataset, batch_size=batch_size, shuffle=True, 
                                   num_workers=num_workers, pin_memory=False)
        hw_val_loader = DataLoader(hw_val_dataset, batch_size=batch_size, shuffle=False, 
                                 num_workers=num_workers, pin_memory=False)
        hw_test_loader = DataLoader(hw_test_dataset, batch_size=batch_size, shuffle=False, 
                                  num_workers=num_workers, pin_memory=False)
        
        log_message(f"PathMNIST loaded: Train={len(hw_train_dataset)}, Val={len(hw_val_dataset)}, Test={len(hw_test_dataset)}")
        
        # Save dataset info
        pathmnist_info = {
            "dataset_name": dataset_name,
            "n_classes": n_classes,
            "class_names": class_names,
            "train_samples": len(hw_train_dataset),
            "val_samples": len(hw_val_dataset),
            "test_samples": len(hw_test_dataset),
            "batch_size": batch_size,
            "device": DEVICE
        }
        save_json(pathmnist_info, "pathmnist_info.json", "predictions")
        
        # Visualize samples
        hw_samples_viz, hw_labels_viz = next(iter(hw_train_loader))
        fig, axes = plt.subplots(1, 5, figsize=(15, 3))
        for i in range(5):
            img = hw_samples_viz[i].permute(1, 2, 0).cpu().numpy()
            # Unnormalize for visualization
            mean = np.array([0.485, 0.456, 0.406])
            std = np.array([0.229, 0.224, 0.225])
            img = std * img + mean
            img = np.clip(img, 0, 1)
            axes[i].imshow(img)
            axes[i].set_title(f"Label: {class_names[hw_labels_viz[i].squeeze().item()]}")
            axes[i].axis('off')
        plt.suptitle("Sample Images from PathMNIST Training Set")
        plt.tight_layout()
        save_figure(fig, "pathmnist_samples.png", "visualizations")
        
        log_message("PathMNIST homework visualization saved")
        
        return {
            "train_loader": hw_train_loader,
            "val_loader": hw_val_loader,
            "test_loader": hw_test_loader,
            "class_names": class_names,
            "n_classes": n_classes
        }
        
    except Exception as e:
        log_message(f"Error in PathMNIST homework: {e}")
        return None

def main():
    """Main execution function with comprehensive error handling"""
    log_message("=== Starting Lung Nodules XAI Analysis ===")
    
    try:
        # Setup data
        ref_loader, inf_loader, normalize, class_names, logit2name = setup_lung_nodule_data()
        
        # Load model
        model = load_pretrained_model()
        
        # Get sample data for analysis
        samples, labels = next(iter(ref_loader))
        
        # Visualize samples
        visualize_samples_and_save(samples, labels, class_names, "nodule_samples")
        
        # Make predictions
        proc_samples, predictions, probabilities, predicted_classes = make_predictions_and_save(
            model, samples, labels, class_names, logit2name
        )
        
        # Setup XAI tools
        cdam_tool, attention_tool = setup_xai_tools(model)
        
        # Generate explanations
        if attention_tool:
            attention_maps = generate_attention_maps(
                attention_tool, proc_samples, class_names, 
                [class_names[labels[i].item()] for i in range(len(labels))],
                predicted_classes
            )
        
        if cdam_tool:
            target_indices = [0] * len(proc_samples)  # Target class indices
            cdam_maps = generate_cdam_maps(cdam_tool, proc_samples, target_indices)
        
        # Setup outlier detection
        outlier_detector = setup_outlier_detection(ref_loader)
        
        # Run homework
        homework_results = run_homework_pathmnist()
        
        # Create summary report
        summary = {
            "execution_summary": {
                "timestamp": datetime.now().isoformat(),
                "device": DEVICE,
                "nodule_prediction_accuracy": sum(p == t for p, t in zip(predicted_classes, [class_names[labels[i].item()] for i in range(len(labels))])) / len(predicted_classes),
                "attention_maps_generated": attention_tool is not None,
                "cdam_maps_generated": cdam_tool is not None,
                "outlier_detection_setup": outlier_detector is not None,
                "homework_completed": homework_results is not None
            },
            "files_generated": {
                "predictions": "predictions/predictions_results.json",
                "dataset_info": "dataset_info.json",
                "attention_maps": "explanations/attention_maps.npy" if attention_tool else None,
                "cdam_maps": "explanations/cdam_maps.npy" if cdam_tool else None,
                "pathmnist_info": "predictions/pathmnist_info.json" if homework_results else None
            }
        }
        
        save_json(summary, "execution_summary.json")
        
        log_message("=== Analysis Complete ===")
        log_message(f"All results saved to: {RESULTS_DIR.absolute()}")
        
        return summary
        
    except Exception as e:
        log_message(f"Critical error in main execution: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Check dependencies
    log_message("Checking dependencies...")
    
    try:
        import obzai
        log_message("obzai package available")
    except ImportError:
        log_message("obzai package not found. Install with: pip install obzai")
        print("Please install required packages:")
        print("pip install torch torchvision transformers medmnist gdown matplotlib numpy obzai")
        sys.exit(1)
    
    # Run main analysis
    results = main()
    
    if results:
        print(f"\n✅ Analysis completed successfully!")
        print(f"📁 Results saved to: {RESULTS_DIR.absolute()}")
        print(f"📊 Check execution_summary.json for overview")
    else:
        print("❌ Analysis failed. Check logs for details.")
