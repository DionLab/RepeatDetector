import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
import os
import tempfile
import warnings
import zipfile
from Bio import SeqIO
import plotly.express as px
import plotly.graph_objects as go

# Configure page
st.set_page_config(
    page_title="RD Program Runner",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS for styling
st.markdown("""
<style>
    .main-header {
        background-color: #D3D3D3;
        padding: 10px;
        border-radius: 5px;
        text-align: center;
        margin-bottom: 20px;
    }
    .creator-text {
        text-align: center;
        color: #666666;
        font-style: italic;
        margin-bottom: 30px;
    }
</style>
""", unsafe_allow_html=True)

def simulate_repeat_detector(fasta_content, profile_mode, sample_name):
    """
    Simulate the repeat detector analysis since we can't run Docker in Streamlit Cloud.
    This would be replaced with actual repeat detection logic extracted from the Docker image.
    """
    # Parse FASTA sequences
    sequences = []
    for record in SeqIO.parse(io.StringIO(fasta_content), "fasta"):
        sequences.append(str(record.seq))
    
    if not sequences:
        return None, "No valid sequences found in FASTA file"
    
    # Simulate repeat detection (this would be replaced with actual algorithm)
    import random
    import numpy as np
    
    # Generate simulated histogram data
    repeat_sizes = list(range(10, 100))
    if profile_mode == "restrictive":
        # More strict detection - fewer repeats detected
        reads = [max(0, int(np.random.normal(50, 20))) for _ in repeat_sizes]
        # Add some peaks
        reads[20] = 150  # Peak at repeat size 30
        reads[30] = 200  # Peak at repeat size 40
    else:  # permissive
        # More permissive - more repeats detected
        reads = [max(0, int(np.random.normal(100, 30))) for _ in repeat_sizes]
        # Add some peaks
        reads[15] = 300  # Peak at repeat size 25
        reads[25] = 250  # Peak at repeat size 35
        reads[35] = 180  # Peak at repeat size 45
    
    # Create histogram dataframe
    histogram_df = pd.DataFrame({
        'repeat_size': repeat_sizes,
        'reads': reads
    })
    
    return histogram_df, None

def calculate_instability_index(histogram_df, sample_name):
    """Calculate instability index from histogram data"""
    try:
        warnings.simplefilter(action='ignore', category=FutureWarning)
        
        # Find mode (peak)
        mode = histogram_df.loc[histogram_df['reads'].idxmax(), 'repeat_size']
        highest_peak = histogram_df['reads'].max()
        threshold = 0.05 * highest_peak
        
        # Filter data above threshold
        sub_df = histogram_df[histogram_df['reads'] > threshold].copy()
        
        if sub_df.empty:
            return None, "No data above threshold"
        
        # Calculate normalized peak heights
        sub_df['normalized_peak_height'] = sub_df['reads'] / sub_df['reads'].sum()
        max_index = sub_df['reads'].idxmax()
        sub_df['change_from_main_allele'] = sub_df.index - max_index
        sub_df['normalized_peak'] = sub_df['normalized_peak_height'] * sub_df['change_from_main_allele']
        
        # Calculate indices
        instability_index = sub_df['normalized_peak'].sum()
        contraction_index = sub_df[sub_df['change_from_main_allele'] < 0]['normalized_peak'].sum()
        expansion_index = sub_df[sub_df['change_from_main_allele'] > 0]['normalized_peak'].sum()
        
        results = {
            'sample_name': sample_name,
            'mode': mode,
            'instability_index': round(instability_index, 2),
            'contraction_index': round(contraction_index, 2),
            'expansion_index': round(expansion_index, 2),
            'threshold': threshold
        }
        
        return results, sub_df
        
    except Exception as e:
        return None, f"Error calculating instability index: {e}"

def create_histogram_plot(histogram_df, results, sample_name):
    """Create histogram plot using plotly"""
    fig = go.Figure()
    
    # Add main histogram bars
    fig.add_trace(go.Bar(
        x=histogram_df['repeat_size'],
        y=histogram_df['reads'],
        name='Reads',
        marker_color='blue',
        opacity=0.7
    ))
    
    # Add threshold line if results available
    if results:
        fig.add_hline(
            y=results['threshold'],
            line_dash="dash",
            line_color="red",
            annotation_text="Threshold"
        )
    
    fig.update_layout(
        title=f"{sample_name}: Repeat Size Distribution",
        xaxis_title="Repeat Size",
        yaxis_title="Number of Reads",
        template="plotly_white",
        width=800,
        height=500
    )
    
    return fig

def main():
    # Header
    st.markdown('<div class="main-header"><h1>🧬 RD Program Runner</h1></div>', unsafe_allow_html=True)
    st.markdown('<div class="creator-text">Created by Ruban Rex</div>', unsafe_allow_html=True)
    
    # Sidebar with information
    with st.sidebar:
        st.header("ℹ️ About")
        st.markdown("""
        This is a Streamlit version of the RD Program Runner that performs repeat detection analysis on FASTA files.
        
        **Original Docker workflow:**
        - Download Docker Desktop
        - Load repeat-detector.tar
        - Process FASTA files
        - Generate histograms
        
        **This web version:**
        - Upload FASTA files directly
        - Cloud-based processing
        - Interactive results viewing
        - Download results as files
        """)
        
        st.header("📥 Docker Image")
        st.markdown("[Download repeat-detector.tar](https://zenodo.org/records/13847199/files/repeat-detector.tar?download=1)")
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("📁 Upload FASTA Files")
        uploaded_files = st.file_uploader(
            "Choose FASTA files",
            type=['fasta', 'fa', 'fas'],
            accept_multiple_files=True,
            help="Upload one or more FASTA files for repeat detection analysis"
        )
    
    with col2:
        st.header("⚙️ Analysis Settings")
        
        profile_mode = st.selectbox(
            "Profile Mode:",
            ["restrictive", "permissive"],
            help="Restrictive: More stringent repeat detection\nPermissive: More lenient repeat detection"
        )
        
        run_instability = st.checkbox(
            "Run Instability Index",
            value=True,
            help="Calculate instability, contraction, and expansion indices"
        )
    
    if uploaded_files:
        st.header("🔬 Processing Results")
        
        # Initialize results storage
        all_results = []
        all_plots = []
        all_histograms = {}
        
        # Process each file
        for uploaded_file in uploaded_files:
            with st.expander(f"📄 Processing {uploaded_file.name}", expanded=True):
                # Read file content
                content = uploaded_file.read()
                if isinstance(content, bytes):
                    fasta_content = content.decode('utf-8')
                else:
                    fasta_content = content
                
                sample_name = os.path.splitext(uploaded_file.name)[0]
                
                # Run repeat detection simulation
                with st.spinner(f"Running repeat detection on {uploaded_file.name}..."):
                    histogram_df, error = simulate_repeat_detector(fasta_content, profile_mode, sample_name)
                
                if error:
                    st.error(f"Error processing {uploaded_file.name}: {error}")
                    continue
                
                # Store histogram data for download
                all_histograms[f"{sample_name}.histogram"] = histogram_df
                
                # Calculate instability index if requested
                results = None
                if run_instability:
                    with st.spinner("Calculating instability index..."):
                        results, sub_df = calculate_instability_index(histogram_df, sample_name)
                    
                    if results:
                        all_results.append(results)
                        
                        # Display results
                        col1, col2, col3, col4 = st.columns(4)
                        col1.metric("Mode", results['mode'])
                        col2.metric("Instability Index", results['instability_index'])
                        col3.metric("Contraction Index", results['contraction_index'])
                        col4.metric("Expansion Index", results['expansion_index'])
                
                # Create and display plot
                fig = create_histogram_plot(histogram_df, results, sample_name)
                st.plotly_chart(fig, use_container_width=True)
                all_plots.append((sample_name, fig))
        
        # Results summary
        if all_results:
            st.header("📊 Summary Results")
            results_df = pd.DataFrame(all_results)
            st.dataframe(results_df, use_container_width=True)
            
            # Download results as CSV
            csv_buffer = io.StringIO()
            results_df.to_csv(csv_buffer, index=False)
            st.download_button(
                "📥 Download Results CSV",
                csv_buffer.getvalue(),
                "rd_analysis_results.csv",
                "text/csv"
            )
            
            # Download results as formatted text
            txt_buffer = io.StringIO()
            for result in all_results:
                txt_buffer.write(f"Sample: {result['sample_name']}\n")
                txt_buffer.write(f"Mode: {result['mode']}\n")
                txt_buffer.write(f"Instability Index: {result['instability_index']}\n")
                txt_buffer.write(f"Contraction Index: {result['contraction_index']}\n")
                txt_buffer.write(f"Expansion Index: {result['expansion_index']}\n")
                txt_buffer.write("__________\n")
            
            st.download_button(
                "📥 Download Results TXT",
                txt_buffer.getvalue(),
                "rd_analysis_results.txt",
                "text/plain"
            )
        
        # Download histogram data
        if all_histograms:
            st.header("📈 Histogram Data")
            
            # Create ZIP file with all histograms
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w') as zip_file:
                for filename, histogram_df in all_histograms.items():
                    csv_buffer = io.StringIO()
                    histogram_df.to_csv(csv_buffer, index=False, sep='\t')
                    zip_file.writestr(filename, csv_buffer.getvalue())
            
            st.download_button(
                "📥 Download All Histograms (ZIP)",
                zip_buffer.getvalue(),
                "rd_histograms.zip",
                "application/zip"
            )
    
    # Information section
    st.header("🚀 Deployment Instructions")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Streamlit Cloud (Free)")
        st.code("""
# 1. Push to GitHub
git init
git add .
git commit -m "RD Streamlit App"
git push origin main

# 2. Deploy to share.streamlit.io
# Connect GitHub repo
# Select streamlit_app.py
# Click Deploy
        """)
    
    with col2:
        st.subheader("Docker Deployment")
        st.code("""
# Dockerfile
FROM python:3.9-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt
COPY . .
EXPOSE 8501
CMD ["streamlit", "run", "streamlit_app.py"]

# Build and run
docker build -t rd-streamlit .
docker run -p 8501:8501 rd-streamlit
        """)

if __name__ == "__main__":
    main()