"""
Mueller matrix analysis: constructing Mueller matrices, Hermitian covariance matrix, and computing Indices of Polarimetric Purity (IPPs).
"""
import numpy as np

def assemble_mueller_matrix(outputs):
    """
    Assemble the 4x4 Mueller matrix from output Stokes vectors for a set of known input polarization states.
    The `outputs` dict should map polarization states to Stokes outputs.
    Required keys: 'H', 'V', 'P45', 'R' corresponding to horizontal, vertical, +45°, and right-circular inputs.
    """
    # Define input Stokes basis vectors
    S_H_in = np.array([1, 1, 0, 0], dtype=float)
    S_V_in = np.array([1,-1, 0, 0], dtype=float)
    S_P45_in = np.array([1, 0, 1, 0], dtype=float)
    S_R_in = np.array([1, 0, 0, 1], dtype=float)
    S_in_matrix = np.column_stack((S_H_in, S_V_in, S_P45_in, S_R_in))
    S_out_matrix = np.column_stack((outputs['H'], outputs['V'], outputs['P45'], outputs['R']))
    # Solve for Mueller matrix M such that M * S_in = S_out
    M = np.dot(S_out_matrix, np.linalg.inv(S_in_matrix))
    return M

def hermitian_covariance_matrix(M):
    """
    Compute the Hermitian covariance matrix H from Mueller matrix M (Eq. 7 in paper) [oai_citation:7‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=pure%20elements%20,as%20follows%3A%20H%20%3D%201) [oai_citation:8‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=%E2%8E%9D%20m00%20%2B%20m01%20%2B,%E2%88%92%20m10%20%2B%20m11%20%E2%8E%9E).
    Returns a 4x4 complex Hermitian matrix.
    """
    m = M
    m00, m01, m02, m03 = m[0,0], m[0,1], m[0,2], m[0,3]
    m10, m11, m12, m13 = m[1,0], m[1,1], m[1,2], m[1,3]
    m20, m21, m22, m23 = m[2,0], m[2,1], m[2,2], m[2,3]
    m30, m31, m32, m33 = m[3,0], m[3,1], m[3,2], m[3,3]
    H = np.zeros((4,4), dtype=complex)
    H[0,0] = m00 + m01 + m10 + m11
    H[0,1] = m02 + m12 + 1j*(m03 + m13)
    H[0,2] = m20 + m21 - 1j*(m30 + m31)
    H[0,3] = m22 + m33 + 1j*(m23 - m32)
    H[1,0] = np.conjugate(H[0,1])
    H[1,1] = m00 - m01 + m10 - m11
    H[1,2] = m22 - m33 - 1j*(m23 + m32)
    H[1,3] = m20 - m21 - 1j*(m30 - m31)
    H[2,0] = np.conjugate(H[0,2])
    H[2,1] = np.conjugate(H[1,2])
    H[2,2] = m00 + m01 - m10 - m11
    H[2,3] = m02 - m12 + 1j*(m03 - m13)
    H[3,0] = np.conjugate(H[0,3])
    H[3,1] = np.conjugate(H[1,3])
    H[3,2] = np.conjugate(H[2,3])
    H[3,3] = m00 - m01 - m10 + m11
    H = H * 0.25
    return H

def compute_IPPs(H):
    """
    Compute the Indices of Polarimetric Purity (P1, P2, P3) from the Hermitian matrix H [oai_citation:9‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=The%20eigenvalue%20spectrum%20of%20the,29) [oai_citation:10‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=P1%20%3D%20%CE%BB0%20,%CE%BB0%20%2B%20%CE%BB1%2B%CE%BB2%E2%88%923%CE%BB3%20trH).
    Returns a tuple (P1, P2, P3).
    """
    vals = np.linalg.eigvals(H)
    vals = np.real_if_close(vals)
    vals = np.real(vals)
    vals.sort()
    vals = vals[::-1]  # descending order
    # Ensure 4 eigenvalues
    if len(vals) < 4:
        vals = np.pad(vals, (0, 4-len(vals)), constant_values=0)
    λ0, λ1, λ2, λ3 = vals[0], vals[1], vals[2], vals[3]
    tr = λ0 + λ1 + λ2 + λ3
    if tr == 0:
        return (0.0, 0.0, 0.0)
    P1 = (λ0 - λ1) / tr
    P2 = (λ0 + λ1 - 2*λ2) / tr
    P3 = (λ0 + λ1 + λ2 - 3*λ3) / tr
    return (P1.real, P2.real, P3.real)













